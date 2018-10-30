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
import logging
import rmgpy.constants as constants
from rmgpy.cantherm.common import get_element_mass
from rmgpy.exceptions import InputError
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
        # Open Molpro log file for parsing
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

        return Natoms - 1

    def loadForceConstantMatrix(self):
        """
        Print the force constant matrix by including the print, hessian command in the input file
        """

        F = None

        Natoms = self.getNumberOfAtoms()
        Nrows = Natoms * 3

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read force constant matrix
            if 'Force Constants (Second Derivatives of the Energy) in [a.u.]' in line:
                F = numpy.zeros((Nrows,Nrows), numpy.float64)
                for i in range(int(math.ceil(Nrows / 5.0))):
                    # Header row
                    line = f.readline()
                    # Matrix element rows
                    for j in range(i*5, Nrows):
                        data = f.readline().split()
                        for k in range(len(data)-1):
                            F[j,i*5+k] = float(data[k+1].replace('D', 'E'))
                            F[i*5+k,j] = F[j,i*5+k]
                # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                F *= 4.35974417e-18 / 5.291772108e-11**2
            line = f.readline()
        # Close file when finished
        f.close()

        return F

    def loadGeometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        Molpro .out file. If multiple such geometries are identified, only the
        last is returned.
        """

        symbol, coord, mass, number = [], [], [], []

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Automatically determine the number of atoms
            if 'Current geometry' in line:
                symbol, coord = [], []
                while 'ENERGY' not in line:
                    line = f.readline()
                line = f.readline()
                while line != '\n':
                    data = line.split()
                    symbol.append(str(data[0]))
                    coord.append([float(data[1]), float(data[2]), float(data[3])])
                    line = f.readline()
                line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()

        #If no optimized coordinates were found, uses the input geometry (for example if reading the geometry from a frequency file
        if coord == []:
            f = open(self.path, 'r')
            line = f.readline()
            while line != '':
                if 'Atomic Coordinates' in line:
                    symbol = []; coord = []
                    for i in range(4):
                        line = f.readline()
                    while line != '\n':
                        data = line.split()
                        symbol.append(str(data[1]))
                        coord.append([float(data[3]), float(data[4]), float(data[5])])
                        line = f.readline()
                line = f.readline()

        # Assign appropriate mass to each atom in the molecule
        for atom1 in symbol:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            number.append(num1)
        number = numpy.array(number, numpy.int)
        mass = numpy.array(mass, numpy.float64)
        coord = numpy.array(coord, numpy.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise InputError('Unable to read atoms from Molpro geometry output file {0}'.format(self.path))

        return coord, number, mass

    def loadConformer(self, symmetry=None, spinMultiplicity=0, opticalIsomers=1, symfromlog=None, label=''):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a MolPro "Freq" quantum chemistry calculation with the thermo printed.
        """

        modes = []
        E0 = 0.0

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':

            # Read the spin multiplicity if not explicitly given
            if spinMultiplicity == 0 and 'spin' in line:
                splits = line.replace('=', ' ').replace(',', ' ').split(' ')
                for i, s in enumerate(splits):
                    if 'spin' in s:
                        spinMultiplicity = int(splits[i+1]) + 1
                        logging.debug(
                            'Conformer {0} is assigned a spin multiplicity of {1}'.format(label, spinMultiplicity))
                        break
            if spinMultiplicity == 0 and 'SPIN SYMMETRY' in line:
                spin_symmetry = line.split()[-1]
                if spin_symmetry == 'Singlet':
                    spinMultiplicity = 1
                elif spin_symmetry == 'Doublet':
                    spinMultiplicity = 2
                elif spin_symmetry == 'Triplet':
                    spinMultiplicity = 3
                elif spin_symmetry == 'Quartet':
                    spinMultiplicity = 4
                elif spin_symmetry == 'Quintet':
                    spinMultiplicity = 5
                elif spin_symmetry == 'Sextet':
                    spinMultiplicity = 6
                if spinMultiplicity:
                    logging.debug(
                        'Conformer {0} is assigned a spin multiplicity of {1}'.format(label, spinMultiplicity))
                    break

            # The data we want is in the Thermochemistry section of the output
            if 'THERMODYNAMICAL' in line:
                modes = []
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
                        if symfromlog is True:
                            symmetry = int(float(line.split()[3]))

                    # Read moments of inertia for external rotational modes
                    elif 'Rotational Constants' in line and line.split()[-1]=='[GHz]':
                        inertia = [float(d) for d in line.split()[-4:-1]]
                        for i in range(3):
                            inertia[i] = constants.h / (8 * constants.pi * constants.pi * inertia[i] * 1e9) *constants.Na*1e23
                        rotation = NonlinearRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                        modes.append(rotation)

                    elif 'Rotational Constant' in line and line.split()[3]=='[GHz]':
                        inertia = float(line.split()[2])
                        inertia = constants.h / (8 * constants.pi * constants.pi * inertia * 1e9) *constants.Na*1e23
                        rotation = LinearRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
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
                if 'vtz' in line.lower() or 'vdz' in line.lower():
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
        logging.debug('Molpro energy found is {0} hartree'.format(E0))
        #multiply E0 by correct constants
        if E0 is not None:
            E0 = E0 * constants.E_h * constants.Na
            logging.debug('Molpro energy found is {0} J/mol'.format(E0))
            return E0
        else: raise Exception('Unable to find energy in Molpro log file.')

    def loadZeroPointEnergy(self):
        """
        Load the unscaled zero-point energy in J/mol from a MolPro log file.
        """

        ZPE = None

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
            raise Exception('Unable to find zero-point energy in MolPro log file. Make sure that the keyword {frequencies, thermo, print,thermo} is included in the input file')

    def loadNegativeFrequency(self):
        """
        Return the negative frequency from a transition state frequency calculation in cm^-1.
        """
        frequency = None
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read vibrational frequencies
            if 'Normal Modes of imaginary frequencies' in line:
                for i in range(3):
                    line = f.readline()
                frequency = line.split()[2]
            line = f.readline()
        f.close()
        if frequency is None:
            raise Exception('Unable to find imaginary frequency in Molpro output file {0}'.format(self.path))
        negativefrequency = -float(frequency)
        return negativefrequency
