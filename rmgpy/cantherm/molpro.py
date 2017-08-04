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

import math
import numpy
import os.path
import rmgpy.constants as constants

from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

class MolproLog:
    """
    Represents a Molpro log file. The attribute `path` refers to the
    location on disk of the Molpro log file of interest. Methods are provided
    to extract a variety of information into CanTherm classes and/or NumPy
    arrays. 
    """
    
    def __init__(self, path):
        self.path = path

    def getNumberOfAtoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the MolPro log file.
        """

        Natoms = 0
        # Open Gaussian log file for parsing
        f = open(self.path, 'r')
        line = f.readline()
        while line != '' and Natoms == 0:
            # Automatically determine the number of atoms
            if 'ATOMIC COORDINATES' in line and Natoms == 0:
                for i in range(4): line = f.readline()
                while 'Bond lengths' not in line:
                    Natoms += 1
                    line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()
        # Return the result
        Natoms -= 1

        return Natoms

    def loadForceConstantMatrix(self):
        """
        No force constant matrices are reported by the MolPro Files
        """

        F = None

        return F

    def loadGeometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Gaussian log file. If multiple such geometries are identified, only the
        last is returned.
        """

        symbol = []; coord = []

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Automatically determine the number of atoms
            if 'Current geometry' in line:
                number = []; coord = []
                for i in range(4): line = f.readline()
                count = 0
                while line != '\n':
                    data = line.split()
                    symbol.append(str(data[0]))
                    coord.append([float(data[1]), float(data[2]), float(data[3])])
                    count += 1
                    line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()

        coord = numpy.array(coord, numpy.float64)
        number = numpy.zeros(len(symbol), numpy.int)
        mass = numpy.zeros(len(symbol), numpy.float64)
        # Use the atomic mass of the most common isotope rather than the
        # average atomic mass
        # These values were taken from "Atomic Weights and Isotopic Compositions" v3.0 (July 2010) from NIST
        for i in range(len(symbol)):
            if symbol[i] == 'H':
                number[i] = 1
                mass[i] = 1.00782503207
            elif symbol[i] == 'C':
                number[i] = 6
                mass[i] = 12.0
            elif symbol[i] == 'N':
                number[i] = 7
                mass[i] = 14.0030740048
            elif symbol[i] == 'O':
                number[i] = 8
                mass[i] = 15.99491461956
            elif symbol[i] == 'P':
                number[i] = 15
                mass[i] = 30.97376163
            elif symbol[i] == 'S':
                number[i] = 16
                mass[i] = 31.97207100
            elif symbol[i] == 'Cl':
                number[i] = 17
                mass[i] = 35.4527
            elif symbol[i] == 'I':
                number[i] = 53
                mass[i] = 126.90447
            else:
                print 'Atomic number {0:d} not yet supported in loadGeometry().'.format(number[i])
        return coord, number, mass

    def loadConformer(self, symmetry=None, spinMultiplicity=None, opticalIsomers=1):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a MolPro "Freq" quantum chemistry calculation with the thermo printed.
        """

        modes = []
        E0 = 0.0

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':

            # The data we want is in the Thermochemistry section of the output
            if 'THERMODYNAMICAL' in line:
                modes = []
                inPartitionFunctions = False
                line = f.readline()
                while line != '':

                    # This marks the end of the thermochemistry section
                    if '*************************************************' in line:
                        break

                    # Read molecular mass for external translational modes
                    elif 'Molecular Mass:' in line:
                        mass = float(line.split()[2])
                        translation = IdealGasTranslation(mass=(mass,"amu"))
                        modes.append(translation)
                    # Read MolPro's estimate of the external symmetry number
                    elif 'Rotational Symmetry factor' in line and symmetry is None:
                        symmetry = int(float(line.split()[3]))

                    # Read moments of inertia for external rotational modes
                    elif 'Rotational Constants' in line and line.split()[-1]=='[GHz]':
                        inertia = [float(d) for d in line.split()[-4:-1]]
                        for i in range(3):
                            inertia[i] = constants.h / (8 * constants.pi * constants.pi * inertia[i] * 1e9) *constants.Na*1e23
                        rotation = NonlinearRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                        modes.append(rotation)

                    elif 'Rotational Constant' in line and line.split()[3]=='[GHz]':
                        inertia = [float(line.split()[2])]
                        inertia[0] = constants.h / (8 * constants.pi * constants.pi * inertia[0] * 1e9) *constants.Na*1e23
                        rotation = LinearRotor(inertia=(inertia[0],"amu*angstrom^2"), symmetry=symmetry)
                        modes.append(rotation)

                    # Read vibrational modes
                    elif 'Vibrational Temperatures' in line:
                        frequencies = []
                        frequencies.extend([float(d) for d in line.split()[3:]])
                        line = f.readline()
                        while line.strip() != '':
                            frequencies.extend([float(d) for d in line.split()])
                            line = f.readline()
                        # Convert from K to cm^-1
                        if len(frequencies) > 0:
                            frequencies = [freq * 0.695039 for freq in frequencies]  # kB = 0.695039 cm^-1/K
                            vibration = HarmonicOscillator(frequencies=(frequencies,"cm^-1"))
                            modes.append(vibration)

                    # Read the next line in the file
                    line = f.readline()

            # Read the next line in the file
            line = f.readline()

        # Close file when finished
        f.close()
        return Conformer(E0=(E0*0.001,"kJ/mol"), modes=modes, spinMultiplicity=spinMultiplicity, opticalIsomers=opticalIsomers)

    def loadEnergy(self,frequencyScaleFactor=1.):
        """
        Return the f12 energy in J/mol from a Molpro Logfile of a CCSD(T)-f12 job.
        This function determines which energy (f12a or f12b) to use based on the basis set,
        which it will parse out of the Molpro file. For the vtz and dtz basis sets f12a is
        better approximation, but for higher basis sets f12b is a better approximation
        """
        f = open(self.path, 'r')
        line=f.readline()
        
        #search for basisSet
        while line!='':
            if 'basis' in line.lower():
                if 'vtz' in line.lower() or'vdz' in line.lower():
                    f12a=True
                else: f12a=False
                break
            line=f.readline()
        else: raise Exception('Could not find basis set in Molpro File')
        #search for energy
        E0=None
        if f12a:
            while line!='':
                if ('RHF-UCCSD(T)-F12a energy' in line
                    or 'RHF-RCCSD(T)-F12a energy' in line
                    or 'CCSD(T)-F12a total energy  ' in line):
                    E0=float(line.split()[-1])
                    break
                if 'Electronic Energy at 0' in line:
                    E0=float(line.split()[-2])
                    break
                line=f.readline()
        else:
            while line!='':
                if ('RHF-UCCSD(T)-F12b energy' in line
                    or 'RHF-RCCSD(T)-F12b energy' in line
                    or 'CCSD(T)-F12b total energy  ' in line):
                    E0=float(line.split()[-1])
                    break
                if 'Electronic Energy at 0' in line:
                    E0=float(line.split()[-2])
                    break
                line=f.readline()

        f.close()
        
        #multiply E0 by correct constants
        if E0 is not None:
            E0 = E0 * constants.E_h * constants.Na
            return E0
        else: raise Exception('Unable to find energy in Molpro log file.')

    def loadZeroPointEnergy(self):
        """
        Load the unscaled zero-point energy in J/mol from a MolPro log file.
        """

        modes = []
        ZPE = None
        spinMultiplicity = 1

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':

            # Do NOT read the ZPE from the "E(ZPE)=" line, as this is the scaled version!
            # We will read in the unscaled ZPE and later multiply the scaling factor
            # from the input file

            if 'Electronic Energy at 0 [K]:' in line:
                electronic_energy = float(line.split()[5])
                line = f.readline()
                EEplusZPE = float(line.split()[5])
                ZPE = (EEplusZPE-electronic_energy) * constants.E_h * constants.Na
            line = f.readline()

        # Close file when finished
        f.close()

        if ZPE is not None:
            return ZPE
        else:
            raise Exception('Unable to find zero-point energy in MolPro log file.')