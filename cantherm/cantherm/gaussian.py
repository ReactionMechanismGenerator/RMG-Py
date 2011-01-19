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
import numpy

from chempy.geometry import Geometry
from chempy.states import Translation, RigidRotor, HinderedRotor, HarmonicOscillator, StatesModel
from chempy import constants

################################################################################

class GaussianLog:
    """
    Represent a log file from Gaussian. The attribute `path` refers to the
    location on disk of the Gaussian log file of interest. Methods are provided
    to extract a variety of information into ChemPy classes and/or NumPy arrays.
    """

    def __init__(self, path):
        self.path = path

    def getNumberOfAtoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the Gaussian log file.
        """

        Natoms = 0
        # Open Gaussian log file for parsing
        f = open(self.path, 'r')
        line = f.readline()
        while line != '' and Natoms == 0:
            # Automatically determine the number of atoms
            if 'Input orientation:' in line and Natoms == 0:
                for i in range(5): line = f.readline()
                while '---------------------------------------------------------------------' not in line:
                    Natoms += 1
                    line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()
        # Return the result
        return Natoms

    def loadForceConstantMatrix(self):
        """
        Return the force constant matrix from the Gaussian log file. The job
        that generated the log file must have the option ``iop(7/33=1)`` in
        order for the proper force constant matrix (in Cartesian coordinates)
        to be printed in the log file. If multiple such matrices are identified,
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
            if 'Force constants in Cartesian coordinates:' in line:
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
        Gaussian log file. If multiple such geometries are identified, only the
        last is returned.
        """

        number = []; coord = []

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Automatically determine the number of atoms
            if 'Input orientation:' in line:
                number = []; coord = []
                for i in range(5): line = f.readline()
                while '---------------------------------------------------------------------' not in line:
                    data = line.split()
                    number.append(int(data[1]))
                    coord.append([float(data[3]), float(data[4]), float(data[5])])
                    line = f.readline()
            line = f.readline()
        # Close file when finished
        f.close()

        coord = numpy.array(coord, numpy.float64)
        number = numpy.array(number, numpy.int)
        mass = numpy.zeros(len(number), numpy.float64)
        # Use the atomic mass of the most common isotope rather than the
        # average atomic mass
        # These values were taken from "Atomic Weights and Isotopic Compositions" v3.0 (July 2010) from NIST
        for i in range(len(number)):
            if number[i] == 1:
                mass[i] = 1.00782503207
            elif number[i] == 6:
                mass[i] = 12.0
            elif number[i] == 8:
                mass[i] = 15.99491461956
            else:
                print 'Atomic number %i not yet supported in loadGeometry().' % number[i]
        
        return Geometry(coordinates=coord * 1e-10, number=number, mass=mass / 1000)

    def loadStates(self, symmetry=None):
        """
        Load the molecular degree of freedom data from a log file created as
        the result of a Gaussian "Freq" quantum chemistry calculation. As
        Gaussian's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value; if
        not provided, the value in the Gaussian log file will be adopted. In a
        log file with multiple Thermochemistry sections, only the last one will
        be kept.
        """

        modes = []
        E0 = 0.0
        spinMultiplicity = 1

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':

            # The data we want is in the Thermochemistry section of the output
            if '- Thermochemistry -' in line:
                modes = []
                inPartitionFunctions = False
                line = f.readline()
                while line != '':

                    # This marks the end of the thermochemistry section
                    if '-------------------------------------------------------------------' in line:
                        break

                    # Read molecular mass for external translational modes
                    elif 'Molecular mass:' in line:
                        mass = float(line.split()[2]) * 1e-3
                        translation = Translation(mass=mass)
                        modes.append(translation)

                    # Read Gaussian's estimate of the external symmetry number
                    elif 'Rotational symmetry number' in line and symmetry is None:
                        symmetry = int(float(line.split()[3]))

                    # Read moments of inertia for external rotational modes
                    elif 'Rotational constants (GHZ):' in line:
                        inertia = [float(d) for d in line.split()[-3:]]
                        for i in range(3):
                            inertia[i] = constants.h / (8 * constants.pi * constants.pi * inertia[i] * 1e9)
                        rotation = RigidRotor(linear=False, inertia=inertia, symmetry=symmetry)
                        modes.append(rotation)
                    elif 'Rotational constant (GHZ):' in line:
                        inertia = [float(line.split()[3])]
                        inertia[0] = constants.h / (8 * constants.pi * constants.pi * inertia[0] * 1e9)
                        rotation = RigidRotor(linear=True, inertia=inertia, symmetry=symmetry)
                        modes.append(rotation)

                    # Read vibrational modes
                    elif 'Vibrational temperatures:' in line:
                        frequencies = []
                        frequencies.extend([float(d) for d in line.split()[2:]])
                        line = f.readline()
                        frequencies.extend([float(d) for d in line.split()[1:]])
                        line = f.readline()
                        while line.strip() != '':
                            frequencies.extend([float(d) for d in line.split()])
                            line = f.readline()
                        # Convert from K to cm^-1
                        frequencies = [freq * 0.695039 for freq in frequencies]  # kB = 0.695039 cm^-1/K
                        vibration = HarmonicOscillator(frequencies=frequencies)
                        modes.append(vibration)

                    # Read ground-state energy
                    elif 'Sum of electronic and zero-point Energies=' in line:
                        E0 = float(line.split()[6]) * 4.35974394e-18 * constants.Na

                    # Read spin multiplicity
                    elif 'Electronic' in line and inPartitionFunctions:
                        spinMultiplicity = int(float(line.split()[1].replace('D', 'E')))

                    elif 'Log10(Q)' in line:
                        inPartitionFunctions = True

                    # Read the next line in the file
                    line = f.readline()

            # Read the next line in the file
            line = f.readline()

        # Close file when finished
        f.close()

        return StatesModel(modes=modes, spinMultiplicity=spinMultiplicity)

    def loadEnergy(self):
        """
        Load the energy in J/mol from a Gaussian log file. The file is checked 
        for a complete basis set extrapolation; if found, that value is 
        returned. Only the last energy in the file is returned.
        """

        modes = []
        E0 = None; E0_cbs = None
        spinMultiplicity = 1

        f = open(self.path, 'r')
        line = f.readline()
        while line != '':

            if 'SCF Done:' in line:
                E0 = float(line.split()[4]) * 4.35974394e-18 * constants.Na
            elif 'CBS-QB3 (0 K)' in line or 'G3 (O K)' in line:
                E0_cbs = float(line.split()[3]) * 4.35974394e-18 * constants.Na
            
            # Read the next line in the file
            line = f.readline()

        # Close file when finished
        f.close()

        if E0_cbs is not None: return E0_cbs
        elif E0 is not None: return E0
        else: raise ChemPyError('Unable to find energy in Gaussian log file.')
    
    def loadScanEnergies(self):
        """
        Extract the optimized energies in J/mol from a log file, e.g. the result
        of a Gaussian "Scan" quantum chemistry calculation.
        """

        # The array of potentials at each scan angle
        Vlist = []

        # Parse the Gaussian log file, extracting the energies of each
        # optimized conformer in the scan
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # The lines containing "SCF Done" give the energy at each
            # iteration (even the intermediate ones)
            if 'SCF Done:' in line:
                E = float(line.split()[4])
            # We want to keep the values of E that come most recently before
            # the line containing "Optimization completed", since it refers
            # to the optimized geometry
            if 'Optimization completed.' in line:
                Vlist.append(E)
            line = f.readline()
        # Close file when finished
        f.close()

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        Vlist = numpy.array(Vlist, numpy.float64)
        Vlist -= numpy.min(Vlist)
        Vlist *= 4.35974394e-18 * 6.02214179e23

        return Vlist

    def loadNegativeFrequency(self):
        """
        Return the negative frequency from a transition state frequency
        calculation in cm^-1.
        """
        
        frequencies = []
        
        f = open(self.path, 'r')
        line = f.readline()
        while line != '':
            # Read vibrational frequencies
            if 'Frequencies --' in line:
                frequencies.extend(line.split()[2:])
            line = f.readline()
        # Close file when finished
        f.close()
        
        frequencies = [float(freq) for freq in frequencies]
        frequencies.sort()
        frequency = [freq for freq in frequencies if freq < 0][0]
        
        return frequency

    def fitCosinePotential(self):
        """
        For a given log file, extract the energies and fit them to a potential
        of the form

        .. math:: V(\\phi) = \\frac{1}{2} V_0 \\left( 1 - \\cos \\sigma \\phi \\right)

        This function returns the fitted barrier height :math:`V_0` in J/mol  
        and symmetry number (degeneracy) :math:`\\sigma`.
        For best results, the Scan should only be performed on one rotor at a
        time. It should begin from the minimum energy conformation, cover one
        complete rotation, and return to the minimum energy conformation as the
        last step in the scan.
        """
        
        # Load the energies from the file
        Vlist = self.loadScanEnergies()
        
        # Gaussian does something extra with the last step in the scan, so we
        # discard this point
        Vlist = Vlist[:-1]

        # Determine the set of dihedral angles corresponding to the above
        # This assumes that you start at 0.0, finish at 360.0, and take
        # constant step sizes in between
        angle = numpy.arange(0.0, 2*math.pi+0.00001, 2*math.pi/(len(Vlist)-1), numpy.float64)

        # Fit the simple cosine potential to get the barrier height V0
        # and the symmetry number
        # We fit at integral symmetry numbers in the range [1, 9]
        # The best fit will have the maximum barrier height
        symmetry = 0; barrier = 0.0
        for symm in range(1, 10):
            num = numpy.sum(Vlist * (1 - numpy.cos(symm * angle)))
            den = numpy.sum((1 - numpy.cos(symm * angle))**2)
            V = 2 * num / den
            if V > barrier:
                symmetry = symm
                barrier = V

        return barrier, symmetry

    def fitFourierSeriesPotential(self):
        """
        For a given log file, extract the energies and fit them to a potential
        of the form

        .. math:: V(\\phi) = \\sum_{m=1}^5 A_m \\cos m \\phi + \\sum_{m=1}^5 B_m \\sin m \\phi
        
        This function returns the fitted Fourier coefficients :math:`A_m` and
        :math:`B_m` in J/mol.
        For best results, the Scan should only be performed on one rotor at a
        time. It should begin from the minimum energy conformation, cover one
        complete rotation, and return to the minimum energy conformation as the
        last step in the scan.
        """

        # Load the energies from the file
        Vlist = self.loadScanEnergies()

        # Gaussian does something extra with the last step in the scan, so we
        # discard this point
        Vlist = Vlist[:-1]
        
        # Determine the set of dihedral angles corresponding to the above
        # This assumes that you start at 0.0, finish at 360.0, and take
        # constant step sizes in between
        angle = numpy.arange(0.0, 2*math.pi+0.00001, 2*math.pi/(len(Vlist)-1), numpy.float64)

        # Fit Fourier series potential
        A = numpy.zeros((len(Vlist),12), numpy.float64)
        b = numpy.zeros(len(Vlist), numpy.float64)
        for i in range(len(Vlist)):
            for m in range(6):
                A[i,m] = math.cos(m * angle[i])
                A[i,6+m] = math.sin(m * angle[i])
                b[i] = Vlist[i]
        x, residues, rank, s = numpy.linalg.lstsq(A, b)

        # Return the set of Fourier coefficients
        return numpy.array([x[1:6], x[7:]], numpy.float64)

################################################################################
