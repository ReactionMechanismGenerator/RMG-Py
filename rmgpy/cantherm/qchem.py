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
import os.path
import rmgpy.constants as constants
from rmgpy.cantherm.common import checkConformerEnergy
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer
################################################################################
class QchemLog:
    """
    Represent an output file from Qchem. The attribute `path` refers to the
    location on disk of the Qchem output file of interest. Methods are provided
    to extract a variety of information into CanTherm classes and/or NumPy
    arrays.
    """    

    def __init__(self, path):
        self.path = path

    def getNumberOfAtoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Qchem output file.
        """

        Natoms = 0
        # Open Qchem log file for parsing
        f = open(self.path, 'r')
        line = f.readline()
        while line != '' and Natoms == 0:
            # Automatically determine the number of atoms
            if 'Standard Nuclear Orientation' in line and Natoms == 0:
                for i in range(3): line = f.readline()
                while '----------------------------------------------------' not in line:
                    Natoms += 1
                    line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()
        # Return the result
        return Natoms

    def loadForceConstantMatrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the 
        QChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """

        F = None

        Natoms = self.getNumberOfAtoms()
        Nrows = Natoms * 3
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read force constant matrix
            if 'Final Hessian.' in line or 'Hessian of the SCF Energy' in line:
                F = numpy.zeros((Nrows,Nrows), numpy.float64)
                for i in range(int(math.ceil(Nrows / 6.0))):
                    # Header row
                    line = f.readline()
                    # Matrix element rows
                    for j in range(Nrows): #for j in range(i*6, Nrows):
                        data = f.readline().split()
                        for k in range(len(data)-1):
                            F[j,i*6+k] = float(data[k+1])
                            #F[i*5+k,j] = F[j,i*5+k]
                # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                F *= 4.35974417e-18 / 5.291772108e-11**2
            line = f.readline()
        # Close file when finished
        f.close()

        return F        
    
    def loadGeometry(self):
        
        """
        Return the optimum geometry of the molecular configuration from the
        Qchem log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atom = []; coord = []; number = []; 

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            if 'Final energy is' in line:
                print 'found a sucessfully completed Qchem Geometry Optimization Job'
                line = f.readline()
                atom = []; coord = []
                break
            line = f.readline()
        found = 0           
        while line != '':        
            if 'Standard Nuclear Orientation' in line:
                found += 1
                for i in range(3): line = f.readline() # skip  lines
                while '----------------------------------------------------' not in line:
                    data = line.split()
                    atom.append((data[1]))
                    coord.append([float(data[2]), float(data[3]), float(data[4])])
                    line = f.readline()
                # Read the next line in the file    
                line = f.readline()
            # Read the next line in the file
            line = f.readline()
            if found ==1: break
        line = f.readline()
        #print coord
        f.close()
        coord = numpy.array(coord, numpy.float64)
        mass = numpy.array(coord, numpy.float64)
        # Assign appropriate mass to each atom in molecule
        # These values were taken from "Atomic Weights and Isotopic Compositions" v3.0 (July 2010) from NIST

        mass = [0]*len(atom)  
        
        for i in range(len(atom)):  
            if atom[i] == 'H':
                mass[i] = 1.00782503207
                number.append('1')
            elif atom[i] == 'C':
                mass[i] = 12.0
                number.append('6')
            elif atom[i] == 'N':
                mass[i] = 14.0030740048
                number.append('7')
            elif atom[i] == 'O':
                mass[i] = 15.99491461956
                number.append('8')
            elif atom[i] == 'P':
                mass[i] = 30.97376163
                number.append('15')
            elif atom[i] == 'S':
                mass[i] = 31.97207100
                number.append('16')
            elif atom[i] == 'Cl':
                mass[i] = 35.4527
                number.append('17')
            else:
                raise NotImplementedError('Atomic atom {0:d} not yet supported in loadGeometry().'.format(atom[i]))
        number = numpy.array(number, numpy.int)       
        return coord, number, mass
    
    def loadConformer(self, symmetry=None, spinMultiplicity=0, opticalIsomers=1, symfromlog=None, label=''):
        """
        Load the molecular degree of freedom data from a output file created as
        the result of a Qchem "Freq"  calculation. As
        Qchem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value; if
        not provided, the value in the Qchem output file will be adopted.
        """
        modes = []; freq = []; mmass = []; rot = []
        E0 = 0.0
#        symmetry = 1
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read spin multiplicity if not explicitly given
            if '$molecule' in line and spinMultiplicity == 0:
                line = f.readline()
                if len(line.split()) == 2:
                    spinMultiplicity = int(float(line.split()[1]))
                    logging.debug('Conformer {0} is assigned a spin multiplicity of {1}'.format(label,spinMultiplicity))
            # The rest of the data we want is in the Thermochemistry section of the output
            elif 'VIBRATIONAL ANALYSIS' in line:
                modes = []

                line = f.readline()
                while line != '':

                    # This marks the end of the thermochemistry section
                    if 'Thank you very much for using Q-Chem.' in line:
                        break

                    # Read vibrational modes
                    elif 'VIBRATIONAL FREQUENCIES (CM**-1)' in line:
                        frequencies = []
                        while 'STANDARD THERMODYNAMIC QUANTITIES AT' not in line:
                            if ' Frequency:' in line:
                                if len(line.split()) == 4:
                                    frequencies.extend([float(d) for d in line.split()[-3:]])
                                elif len(line.split()) == 3:
                                    frequencies.extend([float(d) for d in line.split()[-2:]])
                                elif len(line.split()) == 2:
                                    frequencies.extend([float(d) for d in line.split()[-1:]])    
                            line = f.readline()
                        line = f.readline()
                        # If there is an imaginary frequency, remove it
                        if frequencies[0] < 0.0:
                            frequencies = frequencies[1:]
                            
                        vibration = HarmonicOscillator(frequencies=(frequencies,"cm^-1"))
                        #modes.append(vibration)
                        freq.append(vibration)
                    # Read molecular mass for external translational modes
                    elif 'Molecular Mass:' in line:
                        mass = float(line.split()[2])
                        translation = IdealGasTranslation(mass=(mass,"amu"))
                        #modes.append(translation)
                        mmass.append(translation)

                    # Read moments of inertia for external rotational modes, given in atomic units
                    elif 'Eigenvalues --' in line:
                        inertia = [float(d) for d in line.split()[-3:]]
                        # If the first eigenvalue is 0, the rotor is linear
                        symmetry = 1
                        if inertia[0] == 0.0:
                            inertia.remove(0.0)
                            logging.debug('inertia is {}'.format(str(inertia)))
                            for i in range(2):
                                inertia[i] *= (constants.a0/1e-10)**2
                            inertia = numpy.sqrt(inertia[0]*inertia[1])
                            rotation = LinearRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)    
                            rot.append(rotation)                             
                        else:
                            for i in range(3):
                                inertia[i] *= (constants.a0/1e-10)**2
                                rotation = NonlinearRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                                #modes.append(rotation)
                            rot.append(rotation) 

                    # Read Qchem's estimate of the external rotational symmetry number, which may very well be incorrect
                    elif 'Rotational Symmetry Number is' in line: # and symmetry is None:
                        if symfromlog is True:
                            symmetry = int(float(line.split()[4]))
                            logging.debug('rot sym is {}'.format(str(symmetry)))

                    # Read the next line in the file
                    line = f.readline()

            # Read the next line in the file
            line = f.readline()

        # Close file when finished
        f.close()
        modes = mmass + rot + freq
        return Conformer(E0=(E0*0.001,"kJ/mol"), modes=modes, spinMultiplicity=spinMultiplicity, opticalIsomers=opticalIsomers)
              
    def loadEnergy(self,frequencyScaleFactor=1.):
        """
        Load the energy in J/mol from a Qchem log file.  Only the last energy 
        in the file is returned. The zero-point energy is *not* included in 
        the returned value.
        """
        E0 = None
    
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
    
            if 'Final energy is' in line:
                E0 = float(line.split()[3]) * constants.E_h * constants.Na
                logging.debug('energy is {}'.format(str(E0)))
            
#            elif 'Zero point vibrational energy' in line:
                #Qchem's ZPE is in kcal/mol
#                ZPE = float(line.split()[4]) * 4184
#                scaledZPE = ZPE * frequencyScaleFactor
#                print 'ZPE is ' + str(ZPE)
            # Read the next line in the file
            line = f.readline()
    
        # Close file when finished
        f.close()

        if E0 is not None:
            return E0
        else:
            raise Exception('Unable to find energy in Qchem output file.')
        
    def loadZeroPointEnergy(self,frequencyScaleFactor=1.):
        """
        Load the unscaled zero-point energy in J/mol from a Qchem output file.
        """

        ZPE = None
    
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
    
#            if 'Final energy is' in line:
#                E0 = float(line.split()[3]) * constants.E_h * constants.Na
#                print 'energy is' + str(E0)
            if 'Zero point vibrational energy' in line:
                #Qchem's ZPE is in kcal/mol
                ZPE = float(line.split()[4]) * 4184
                #scaledZPE = ZPE * frequencyScaleFactor
                logging.debug('ZPE is {}'.format(str(ZPE)))
            # Read the next line in the file
            line = f.readline()
    
        # Close file when finished
        f.close()
        
        if ZPE is not None:
            return ZPE
        else:
            raise Exception('Unable to find zero-point energy in Qchem output file.')
              
    def loadScanEnergies(self):
        """
        Extract the optimized energies in J/mol from a Qchem log file, e.g. the 
        result of a Qchem "PES Scan" quantum chemistry calculation.
        """
        Vlist = []
        angle = []
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            if 'Summary of potential scan:' in line:
                line = f.readline()
                print 'found a sucessfully completed Qchem Job'
                while '-----------------' not in line:
        #            print len(line.split())
        #            Vlist.append(float(line.split()[1]))
                    values = [float(item) for item in line.split()]
                    angle.append(values[0])
                    Vlist.append(values[1])
            # Read the next line in the file
                    line = f.readline()
            line = f.readline()
            if 'SCF failed to converge' in line:
                print 'Qchem Job did not sucessfully complete: SCF failed to converge'
                break
        # Close file when finished   
        print '   Assuming', os.path.basename(self.path), 'is the output from a Qchem PES scan...'
        f.close()  
                    

        Vlist = numpy.array(Vlist, numpy.float64)
        # check to see if the scanlog indicates that one of your reacting species may not be the lowest energy conformer
        checkConformerEnergy(Vlist, self.path)
        
        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        Vlist -= numpy.min(Vlist)
        Vlist *= constants.E_h * constants.Na      
        angle = numpy.arange(0.0, 2*math.pi+0.00001, 2*math.pi/(len(Vlist)-1), numpy.float64)
        return Vlist, angle
        
    def loadNegativeFrequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read imaginary frequency
            if ' Frequency:' in line:
                frequency = float((line.split()[1]))
                break
            line = f.readline()
        # Close file when finished
        f.close()
        #Make sure the frequency is imaginary:
        if frequency < 0:
            return frequency
        else:
            raise Exception('Unable to find imaginary frequency in QChem output file {0}'.format(self.path))
