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

import os.path
import numpy
import logging

from rmgpy.cantherm.output import prettify
from rmgpy.kinetics import *

import rmgpy.quantity as quantity
import rmgpy.constants as constants

################################################################################

class KineticsJob:
    """
    A representation of a CanTherm kinetics job. This job is used to compute 
    and save the high-pressure-limit kinetics information for a single reaction.
    """
    
    def __init__(self, reaction,  
                 Tmin=None, 
                 Tmax=None,
                 Tlist=None,
                 Tcount=0):
        
        if Tmin is not None:
            self.Tmin = quantity.Quantity(Tmin)
        else:
            self.Tmin = None
            
        if Tmax is not None:
            self.Tmax = quantity.Quantity(Tmax)
        else:
            self.Tmax = None
        
        self.Tcount = Tcount
        
        if Tlist is not None:
            self.Tlist = quantity.Quantity(Tlist)
            self.Tmin = quantity.Quantity(numpy.min(self.Tlist.value_si),"K")
            self.Tmax = quantity.Quantity(numpy.max(self.Tlist.value_si),"K")
            self.Tcount = len(self.Tlist.value_si)
        else:
            if Tmin and Tmax is not None:
                
                if self.Tcount <= 3.:
                    self.Tcount = 50
                
                stepsize = (self.Tmax.value_si-self.Tmin.value_si)/self.Tcount
                
                self.Tlist = quantity.Quantity(numpy.arange(self.Tmin.value_si, self.Tmax.value_si+stepsize, stepsize),"K")
            else:
                self.Tlist = None
        
        self.reaction = reaction
        self.kunits = None
    
    @property
    def Tmin(self):
        """The minimum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmin
    @Tmin.setter
    def Tmin(self, value):
        self._Tmin = quantity.Temperature(value)
    
    @property
    def Tmax(self):
        """The maximum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmax
    @Tmax.setter
    def Tmax(self, value):
        self._Tmax = quantity.Temperature(value)
    
    @property
    def Tlist(self):
        """The temperatures at which the k(T) values are computed."""
        return self._Tlist
    @Tlist.setter
    def Tlist(self, value):
        self._Tlist = quantity.Temperature(value)
        
    def execute(self, outputFile=None, plot=False):
        """
        Execute the kinetics job, saving the results to the given `outputFile`
        on disk.
        """
        
        if self.Tlist is not None:
            self.generateKinetics(self.Tlist.value_si)
        else:
            self.generateKinetics()
        if outputFile is not None:
            self.save(outputFile)
            if plot:
                self.plot(os.path.dirname(outputFile))
    
    def generateKinetics(self,Tlist=None):
        """
        Generate the kinetics data for the reaction and fit it to a modified
        Arrhenius model.
        """
        kineticsClass = 'Arrhenius'
        
        tunneling = self.reaction.transitionState.tunneling
        if isinstance(tunneling, Wigner) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transitionState.frequency.value_si,"cm^-1")
        elif isinstance(tunneling, Eckart) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transitionState.frequency.value_si,"cm^-1")
            tunneling.E0_reac = (sum([reactant.conformer.E0.value_si for reactant in self.reaction.reactants])*0.001,"kJ/mol")
            tunneling.E0_TS = (self.reaction.transitionState.conformer.E0.value_si*0.001,"kJ/mol")
            tunneling.E0_prod = (sum([product.conformer.E0.value_si for product in self.reaction.products])*0.001,"kJ/mol")
        elif tunneling is not None:
            raise ValueError('Unknown tunneling model {0!r}.'.format(tunneling))
        
        logging.info('Generating {0} kinetics model for {0}...'.format(kineticsClass, self.reaction))
        
        if Tlist is None:
            Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
        klist = numpy.zeros_like(Tlist)
        for i in range(Tlist.shape[0]):
            klist[i] = self.reaction.calculateTSTRateCoefficient(Tlist[i])

        order = len(self.reaction.reactants)
        klist *= 1e6 ** (order-1)
        self.kunits = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
        self.reaction.kinetics = Arrhenius().fitToData(Tlist, klist, kunits=self.kunits)

    def save(self, outputFile):
        """
        Save the results of the kinetics job to the file located
        at `path` on disk.
        """
        reaction = self.reaction
        
        logging.info('Saving kinetics for {0}...'.format(reaction))
        
        order = len(self.reaction.reactants)
        factor = 1e6 ** (order-1)
        
        f = open(outputFile, 'a')
    
        f.write('#   ======= =========== =========== =========== ===============\n')
        f.write('#   Temp.   k (TST)     Tunneling   k (TST+T)   Units\n')
        f.write('#   ======= =========== =========== =========== ===============\n')
        
        if self.Tlist is None:
            Tlist = [300,400,500,600,800,1000,1500,2000]
        else:
            Tlist =self.Tlist.value_si

        for T in Tlist:  
            tunneling = reaction.transitionState.tunneling
            reaction.transitionState.tunneling = None
            k0 = reaction.calculateTSTRateCoefficient(T) * factor
            reaction.transitionState.tunneling = tunneling
            k = reaction.calculateTSTRateCoefficient(T) * factor
            tunneling = reaction.transitionState.tunneling
            kappa = k / k0
            f.write('#    {0:4g} K {1:11.3e} {2:11g} {3:11.3e} {4}\n'.format(T, k0, kappa, k, self.kunits))
        f.write('#   ======= =========== =========== =========== ===============\n')
        
        # Reaction path degeneracy is INCLUDED in the kinetics itself!
        string = 'kinetics(label={0!r}, kinetics={1!r})'.format(reaction.label, reaction.kinetics)
        f.write('{0}\n\n'.format(prettify(string)))
        
        f.close()
        
        # Also save the result to chem.inp
        f = open(os.path.join(os.path.dirname(outputFile), 'chem.inp'), 'a')
        
        reaction = self.reaction
        kinetics = reaction.kinetics
                
        string = '{0!s:51} {1:9.3e} {2:9.3f} {3:9.3f}\n'.format(
            reaction,
            kinetics.A.value_si * factor,
            kinetics.n.value_si,
            kinetics.Ea.value_si / 4184.,
        )
                
        f.write('{0}\n'.format(string))
            
        f.close()
        
    def plot(self, outputDirectory):
        """
        Plot both the raw kinetics data and the Arrhenius fit versus 
        temperature. The plot is saved to the file ``kinetics.pdf`` in the
        output directory. The plot is not generated if ``matplotlib`` is not
        installed.
        """
        # Skip this step if matplotlib is not installed
        try:
            import pylab
        except ImportError:
            return

        Tlist = 1000.0/numpy.arange(0.4, 3.35, 0.05)
        klist = numpy.zeros_like(Tlist)
        klist2 = numpy.zeros_like(Tlist)
        for i in range(Tlist.shape[0]):
            klist[i] = self.reaction.calculateTSTRateCoefficient(Tlist[i])
            klist2[i] = self.reaction.kinetics.getRateCoefficient(Tlist[i])

        order = len(self.reaction.reactants)
        klist *= 1e6 ** (order-1)
        klist2 *= 1e6 ** (order-1)

        pylab.semilogy(1000.0 / Tlist, klist, 'ok')
        pylab.semilogy(1000.0 / Tlist, klist2, '-k')
        pylab.xlabel('1000 / Temperature (1000/K)')
        pylab.ylabel('Rate coefficient ({0})'.format(self.kunits))
        pylab.savefig(os.path.join(outputDirectory, 'kinetics.pdf'))
        pylab.close()
