#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module contains the :class:`ThermoJob` class, used to compute and save the
thermodynamics information for a single species.
"""

import os.path
import numpy.linalg
import logging
import string

import rmgpy.constants as constants
from rmgpy.cantherm.output import prettify

from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor
from rmgpy.statmech.conformer import Conformer

from rmgpy.thermo.thermodata import ThermoData
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.chemkin import writeThermoEntry
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.molecule.util import retrieveElementCount

################################################################################

class ThermoJob(object):
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
        
        if species.thermo is not None:
            logging.info("Thermo already generated for species {}. Skipping thermo generation.".format(species))
            return None
        
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
        H298 = species.getThermoData().getEnthalpy(298) / 4184.
        S298 = species.getThermoData().getEntropy(298) / 4.184
        f.write('#   Enthalpy of formation (298 K)   = {0:9.3f} kcal/mol\n'.format(H298))
        f.write('#   Entropy of formation (298 K)    = {0:9.3f} cal/(mol*K)\n'.format(S298))
        f.write('#    =========== =========== =========== =========== ===========\n')
        f.write('#    Temperature Heat cap.   Enthalpy    Entropy     Free energy\n')
        f.write('#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)\n')
        f.write('#    =========== =========== =========== =========== ===========\n')
        for T in [300,400,500,600,800,1000,1500,2000,2400]:
            try:
                Cp = species.getThermoData().getHeatCapacity(T) / 4.184
                H = species.getThermoData().getEnthalpy(T) / 4184.
                S = species.getThermoData().getEntropy(T) / 4.184
                G = species.getThermoData().getFreeEnergy(T) / 4184.
                f.write('#    {0:11g} {1:11.3f} {2:11.3f} {3:11.3f} {4:11.3f}\n'.format(T, Cp, H, S, G))
            except ValueError:
                logging.debug("Valid thermo for {0} is outside range for temperature {1}".format(species,T))
        f.write('#    =========== =========== =========== =========== ===========\n')
        
        string = 'thermo(label={0!r}, thermo={1!r})'.format(species.label, species.getThermoData())
        f.write('{0}\n\n'.format(prettify(string)))
        
        f.close()
        # write chemkin file
        f = open(os.path.join(os.path.dirname(outputFile), 'chem.inp'), 'a')
        if isinstance(species, Species):
            if species.molecule and isinstance(species.molecule[0], Molecule):
                elementCounts = retrieveElementCount(species.molecule[0])
            else:
                try:
                    elementCounts = species.props['elementCounts']
                except KeyError:
                    elementCounts = {'C': 0, 'H': 0}
        else:
            elementCounts = {'C': 0, 'H': 0}
        string = writeThermoEntry(species, elementCounts=elementCounts, verbose=True)
        f.write('{0}\n'.format(string))
        f.close()

        # write species dictionary
        f = open(os.path.join(os.path.dirname(outputFile), 'species_dictionary.txt'), 'a')
        if isinstance(species, Species):
            if species.molecule and isinstance(species.molecule[0], Molecule):
                f.write(species.molecule[0].toAdjacencyList(removeH=False,label=species.label))
                f.write('\n')
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
            import matplotlib.pyplot as plt
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
        thermo = self.species.getThermoData()
        for i in range(Tlist.shape[0]):
            Cplist[i] = conformer.getHeatCapacity(Tlist[i])
            Slist[i] = conformer.getEntropy(Tlist[i])
            Hlist[i] = (conformer.getEnthalpy(Tlist[i]) + conformer.E0.value_si) * 0.001
            Glist[i] = Hlist[i] - Tlist[i] * Slist[i] * 0.001
            Cplist1[i] = thermo.getHeatCapacity(Tlist[i])
            Slist1[i] = thermo.getEntropy(Tlist[i])
            Hlist1[i] = thermo.getEnthalpy(Tlist[i]) * 0.001
            Glist1[i] = thermo.getFreeEnergy(Tlist[i]) * 0.001

        fig = plt.figure(figsize=(10,8))
        fig.suptitle('{0}'.format(self.species.label))
        plt.subplot(2,2,1)
        plt.plot(Tlist, Cplist / 4.184, '-r', Tlist, Cplist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Heat capacity (cal/mol*K)')
        plt.legend(['statmech', 'fitted'], loc=4)

        plt.subplot(2,2,2)
        plt.plot(Tlist, Slist / 4.184, '-r', Tlist, Slist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Entropy (cal/mol*K)')

        plt.subplot(2,2,3)
        plt.plot(Tlist, Hlist / 4.184, '-r', Tlist, Hlist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Enthalpy (kcal/mol)')

        plt.subplot(2,2,4)
        plt.plot(Tlist, Glist / 4.184, '-r', Tlist, Glist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Gibbs free energy (kcal/mol)')

        fig.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.95, wspace=0.35, hspace=0.20)

        plot_path = os.path.join(outputDirectory, 'plots')

        if not os.path.exists(plot_path):
            os.mkdir(plot_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        filename = ''.join(c for c in self.species.label if c in valid_chars) + '.pdf'
        plt.savefig(os.path.join(plot_path, filename))
        plt.close()
