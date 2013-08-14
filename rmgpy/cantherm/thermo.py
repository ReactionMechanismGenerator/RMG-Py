#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains the :class:`ThermoJob` class, used to compute and save the
thermodynamics information for a single species.
"""

import os.path
import math
import numpy.linalg
import logging

import rmgpy.constants as constants
from rmgpy.cantherm.output import prettify
from rmgpy.statmech import *
from rmgpy.thermo import *

################################################################################

class ThermoJob:
    """
    A representation of a CanTherm thermodynamics job. This job is used to
    compute and save the thermodynamics information for a single species.
    """
    
    def __init__(self, species, thermoClass):
        self.species = species
        self.thermoClass = thermoClass
    
    def execute(self, outputFile=None, plot=False):
        """
        Execute the thermodynamics job, saving the results to the
        given `outputFile` on disk.
        """
        self.generateThermo()
        if outputFile is not None:
            self.save(outputFile)
            if plot:
                self.plot(os.path.dirname(outputFile))
    
    def generateThermo(self):
        """
        Generate the thermodynamic data for the species and fit it to the
        desired heat capacity model (as specified in the `thermoClass` 
        attribute).
        """
        if self.thermoClass.lower() not in ['wilhoit', 'nasa']:
            raise Exception('Unknown thermodynamic model "{0}".'.format(self.thermoClass))
    
        species = self.species
    
        logging.info('Generating {0} thermo model for {1}...'.format(self.thermoClass, species))
        
        Tlist = numpy.arange(10.0, 3001.0, 10.0, numpy.float64)
        Cplist = numpy.zeros_like(Tlist)
        H298 = 0.0
        S298 = 0.0
        conformer = self.species.conformer
        for i in range(Tlist.shape[0]):
            Cplist[i] += conformer.getHeatCapacity(Tlist[i])
        H298 += conformer.getEnthalpy(298.) + conformer.E0.value_si
        S298 += conformer.getEntropy(298.)
        
        if not any([isinstance(mode, (LinearRotor, NonlinearRotor)) for mode in conformer.modes]):
            # Monatomic species
            linear = False
            Nfreq = 0
            Nrotors = 0
            Cp0 = 2.5 * constants.R
            CpInf = 2.5 * constants.R
        else:
            # Polyatomic species
            linear = True if isinstance(conformer.modes[1], LinearRotor) else False
            Nfreq = len(conformer.modes[2].frequencies.value)
            Nrotors = len(conformer.modes[3:])
            Cp0 = (3.5 if linear else 4.0) * constants.R
            CpInf = Cp0 + (Nfreq + 0.5 * Nrotors) * constants.R
    
        wilhoit = Wilhoit()
        if Nfreq == 0 and Nrotors == 0:
            wilhoit.Cp0 = (Cplist[0],"J/(mol*K)") 
            wilhoit.CpInf = (Cplist[0],"J/(mol*K)")
            wilhoit.B = (500.,"K") 
            wilhoit.H0 = (0.0,"J/mol")
            wilhoit.S0 = (0.0,"J/(mol*K)") 
            wilhoit.H0 =  (H298 -wilhoit.getEnthalpy(298.15), "J/mol") 
            wilhoit.S0 = (S298 - wilhoit.getEntropy(298.15),"J/(mol*K)")
        else:
            wilhoit.fitToData(Tlist, Cplist, Cp0, CpInf, H298, S298, B0=500.0)
        
        if self.thermoClass.lower() == 'nasa':
            species.thermo = wilhoit.toNASA(Tmin=10.0, Tmax=3000.0, Tint=500.0)
        else:
            species.thermo = wilhoit

    def save(self, outputFile):
        """
        Save the results of the thermodynamics job to the file located
        at `path` on disk.
        """
        species = self.species
        logging.info('Saving thermo for {0}...'.format(species.label))
        
        f = open(outputFile, 'a')
    
        f.write('# Thermodynamics for {0}:\n'.format(species.label))
        H298 = species.thermo.getEnthalpy(298) / 4184.
        S298 = species.thermo.getEntropy(298) / 4.184
        f.write('#   Enthalpy of formation (298 K)   = {0:9.3f} kcal/mol\n'.format(H298))
        f.write('#   Entropy of formation (298 K)    = {0:9.3f} cal/(mol*K)\n'.format(S298))
        f.write('#    =========== =========== =========== =========== ===========\n')
        f.write('#    Temperature Heat cap.   Enthalpy    Entropy     Free energy\n')
        f.write('#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)\n')
        f.write('#    =========== =========== =========== =========== ===========\n')
        for T in [300,400,500,600,800,1000,1500,2000,2400]:
            Cp = species.thermo.getHeatCapacity(T) / 4.184
            H = species.thermo.getEnthalpy(T) / 4184.
            S = species.thermo.getEntropy(T) / 4.184
            G = species.thermo.getFreeEnergy(T) / 4184.
            f.write('#    {0:11g} {1:11.3f} {2:11.3f} {3:11.3f} {4:11.3f}\n'.format(T, Cp, H, S, G))
        f.write('#    =========== =========== =========== =========== ===========\n')
        
        string = 'thermo(label={0!r}, thermo={1!r})'.format(species.label, species.thermo)
        f.write('{0}\n\n'.format(prettify(string)))
        
        f.close()
        
        f = open(os.path.join(os.path.dirname(outputFile), 'chem.inp'), 'a')
        
        thermo = species.thermo
        if isinstance(thermo, NASA):
        
            poly_low = thermo.polynomials[0]
            poly_high = thermo.polynomials[1]
        
            # Determine the number of each type of element in the molecule
            elements = ['C','H','N','O']; elementCounts = [0,0,0,0]

            # Remove elements with zero count
            index = 2
            while index < len(elementCounts):
                if elementCounts[index] == 0:
                    del elements[index]
                    del elementCounts[index]
                else:
                    index += 1
        
            # Line 1
            string = '{0:<16}        '.format(species.label)
            if len(elements) <= 4:
                # Use the original Chemkin syntax for the element counts
                for symbol, count in zip(elements, elementCounts):
                    string += '{0!s:<2}{1:<3d}'.format(symbol, count)
                string += '     ' * (4 - len(elements))
            else:
                string += '     ' * 4
            string += 'G{0:<10.3f}{1:<10.3f}{2:<8.2f}      1'.format(poly_low.Tmin.value_si, poly_high.Tmax.value_si, poly_low.Tmax.value_si)
            if len(elements) > 4:
                string += '&\n'
                # Use the new-style Chemkin syntax for the element counts
                # This will only be recognized by Chemkin 4 or later
                for symbol, count in zip(elements, elementCounts):
                    string += '{0!s:<2}{1:<3d}'.format(symbol, count)
            string += '\n'
        
            # Line 2
            string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2\n'.format(poly_high.c0, poly_high.c1, poly_high.c2, poly_high.c3, poly_high.c4)
        
            # Line 3
            string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3\n'.format(poly_high.c5, poly_high.c6, poly_low.c0, poly_low.c1, poly_low.c2)
        
            # Line 4
            string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4\n'.format(poly_low.c3, poly_low.c4, poly_low.c5, poly_low.c6)
        
            f.write(string)
            
            f.close()
    

    def plot(self, outputDirectory):
        """
        Plot the heat capacity, enthapy, entropy, and Gibbs free energy of the
        fitted thermodynamics model, along with the same values from the
        statistical mechanics model that the thermodynamics model was fitted 
        to. The plot is saved to the file ``thermo.pdf`` in the output
        directory. The plot is not generated if ``matplotlib`` is not installed.
        """
        # Skip this step if matplotlib is not installed
        try:
            import pylab
        except ImportError:
            return
        
        Tlist = numpy.arange(10.0, 2501.0, 10.0)
        Cplist = numpy.zeros_like(Tlist)
        Cplist1 = numpy.zeros_like(Tlist)
        Hlist = numpy.zeros_like(Tlist)
        Hlist1 = numpy.zeros_like(Tlist)
        Slist = numpy.zeros_like(Tlist)
        Slist1 = numpy.zeros_like(Tlist)
        Glist = numpy.zeros_like(Tlist)
        Glist1 = numpy.zeros_like(Tlist)
        
        conformer = self.species.conformer
        thermo = self.species.thermo
        for i in range(Tlist.shape[0]):
            Cplist[i] = conformer.getHeatCapacity(Tlist[i])
            Slist[i] = conformer.getEntropy(Tlist[i])
            Hlist[i] = (conformer.getEnthalpy(Tlist[i]) + conformer.E0.value_si) * 0.001
            Glist[i] = Hlist[i] - Tlist[i] * Slist[i] * 0.001
            Cplist1[i] = thermo.getHeatCapacity(Tlist[i])
            Slist1[i] = thermo.getEntropy(Tlist[i])
            Hlist1[i] = thermo.getEnthalpy(Tlist[i]) * 0.001
            Glist1[i] = thermo.getFreeEnergy(Tlist[i]) * 0.001

        fig = pylab.figure(figsize=(10,8))

        pylab.subplot(2,2,1)
        pylab.plot(Tlist, Cplist / 4.184, '-r', Tlist, Cplist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Heat capacity (cal/mol*K)')
        pylab.legend(['statmech', 'thermo'], loc=4)

        pylab.subplot(2,2,2)
        pylab.plot(Tlist, Slist / 4.184, '-r', Tlist, Slist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Entropy (cal/mol*K)')

        pylab.subplot(2,2,3)
        pylab.plot(Tlist, Hlist / 4.184, '-r', Tlist, Hlist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Enthalpy (kcal/mol)')

        pylab.subplot(2,2,4)
        pylab.plot(Tlist, Glist / 4.184, '-r', Tlist, Glist1 / 4.184, '-b')
        pylab.xlabel('Temperature (K)')
        pylab.ylabel('Gibbs free energy (kcal/mol)')

        fig.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.95, wspace=0.35, hspace=0.20)
        pylab.savefig(os.path.join(outputDirectory, 'thermo.pdf'))
        pylab.close()
