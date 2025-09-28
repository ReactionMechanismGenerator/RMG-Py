#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Arkane QChem module
Used to parse QChem output files
"""

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
import logging
import os.path
import numpy

import rmgpy.constants as constants
from rmgpy.exceptions import InputError
from rmgpy.statmech import IdealGasTranslation, NonlinearRotor, LinearRotor, HarmonicOscillator, Conformer

from arkane.common import check_conformer_energy, get_element_mass
from arkane.log import Log

################################################################################


class QChemLog(Log):
    """
    Represent an output file from QChem. The attribute `path` refers to the
    location on disk of the QChem output file of interest. Methods are provided
    to extract a variety of information into Arkane classes and/or NumPy
    arrays.
    """

    def __init__(self, path):
        super(QChemLog, self).__init__(path)

    def getNumberOfAtoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the QChem output file.
        """
        Natoms = 0

        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '' and Natoms == 0:
                # Automatically determine the number of atoms
                if 'Standard Nuclear Orientation' in line and Natoms == 0:
                    for i in range(3):
                        line = f.readline()
                    while '----------------------------------------------------' not in line:
                        Natoms += 1
                        line = f.readline()
                line = f.readline()

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
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read force constant matrix
                if 'Final Hessian.' in line or 'Hessian of the SCF Energy' in line:
                    F = numpy.zeros((Nrows, Nrows), numpy.float64)
                    for i in range(int(math.ceil(Nrows / 6.0))):
                        # Header row
                        line = f.readline()
                        # Matrix element rows
                        for j in range(Nrows):  # for j in range(i*6, Nrows):
                            data = f.readline().split()
                            for k in range(len(data) - 1):
                                F[j, i * 6 + k] = float(data[k + 1])
                                # F[i*5+k,j] = F[j,i*5+k]
                    # Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
                    F *= 4.35974417e-18 / 5.291772108e-11 ** 2
                line = f.readline()

        return F

    def loadGeometry(self):

        """
        Return the optimum geometry of the molecular configuration from the
        QChem log file. If multiple such geometries are identified, only the
        last is returned.
        """
        atom, coord, number, mass = [], [], [], []

        with open(self.path) as f:
            log = f.read().splitlines()

        # First check that the QChem job file (not necessarily a geometry optimization)
        # has successfully completed, if not an error is thrown
        completed_job = False
        for line in reversed(log):
            if 'Total job time:' in line:
                logging.debug('Found a sucessfully completed QChem Job')
                completed_job = True
                break

        if not completed_job:
            raise InputError(
                'Could not find a successfully completed QChem job in QChem output file {0}'.format(self.path))

        # Now look for the geometry.
        # Will return the final geometry in the file under Standard Nuclear Orientation.
        geometry_flag = False
        for i in reversed(xrange(len(log))):
            line = log[i]
            if 'Standard Nuclear Orientation' in line:
                for line in log[(i + 3):]:
                    if '------------' not in line:
                        data = line.split()
                        atom.append(data[1])
                        coord.append([float(c) for c in data[2:]])
                        geometry_flag = True
                    else:
                        break
                if geometry_flag:
                    break

        # Assign appropriate mass to each atom in the molecule
        for atom1 in atom:
            mass1, num1 = get_element_mass(atom1)
            mass.append(mass1)
            number.append(num1)
        coord = numpy.array(coord, numpy.float64)
        number = numpy.array(number, numpy.int)
        mass = numpy.array(mass, numpy.float64)
        if len(number) == 0 or len(coord) == 0 or len(mass) == 0:
            raise InputError('Unable to read atoms from QChem geometry output file {0}'.format(self.path))

        return coord, number, mass

    def loadConformer(self, symmetry=None, spinMultiplicity=0, opticalIsomers=None, label=''):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        QChem "Freq" calculation. As QChem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the QChem output file will be adopted.
        """
        modes = []
        freq = []
        mmass = []
        rot = []
        inertia = []
        unscaled_frequencies = []
        e0 = 0.0
        if opticalIsomers is None or symmetry is None:
            _opticalIsomers, _symmetry, _ = self.get_symmetry_properties()
            if opticalIsomers is None:
                opticalIsomers = _opticalIsomers
            if symmetry is None:
                symmetry = _symmetry
        with open(self.path, 'r') as f:
            line = f.readline()
            while line != '':
                # Read spin multiplicity if not explicitly given
                if '$molecule' in line and spinMultiplicity == 0:
                    line = f.readline()
                    if len(line.split()) == 2:
                        spinMultiplicity = int(float(line.split()[1]))
                        logging.debug(
                            'Conformer {0} is assigned a spin multiplicity of {1}'.format(label, spinMultiplicity))
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

                            unscaled_frequencies = frequencies
                            vibration = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
                            # modes.append(vibration)
                            freq.append(vibration)
                        # Read molecular mass for external translational modes
                        elif 'Molecular Mass:' in line:
                            mass = float(line.split()[2])
                            translation = IdealGasTranslation(mass=(mass, "amu"))
                            # modes.append(translation)
                            mmass.append(translation)

                        # Read moments of inertia for external rotational modes, given in atomic units
                        elif 'Eigenvalues --' in line:
                            inertia = [float(d) for d in line.split()[-3:]]

                        # Read the next line in the file
                        line = f.readline()

                # Read the next line in the file
                line = f.readline()

                if len(inertia):
                    if inertia[0] == 0.0:
                        # If the first eigenvalue is 0, the rotor is linear
                        inertia.remove(0.0)
                        logging.debug('inertia is {}'.format(str(inertia)))
                        for i in range(2):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                        inertia = numpy.sqrt(inertia[0] * inertia[1])
                        rotation = LinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                        rot.append(rotation)
                    else:
                        for i in range(3):
                            inertia[i] *= (constants.a0 / 1e-10) ** 2
                            rotation = NonlinearRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                            # modes.append(rotation)
                        rot.append(rotation)

                    inertia = []

        modes = mmass + rot + freq
        return Conformer(E0=(e0 * 0.001, "kJ/mol"), modes=modes, spinMultiplicity=spinMultiplicity,
                         opticalIsomers=opticalIsomers), unscaled_frequencies

    def loadEnergy(self, zpe_scale_factor=1.):
        """
        Load the energy in J/mol from a QChem log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, 'r') as f:
            a = b = 0
            for line in f:
                if 'Final energy is' in line:
                    a = float(line.split()[3]) * constants.E_h * constants.Na
                if 'Total energy in the final basis set' in line:
                    b = float(line.split()[8]) * constants.E_h * constants.Na
                e_elect = a or b
        if e_elect is None:
            raise InputError('Unable to find energy in QChem output file.')
        return e_elect

    def loadZeroPointEnergy(self):
        """
        Load the unscaled zero-point energy in J/mol from a QChem output file.
        """
        zpe = None
        with open(self.path, 'r') as f:
            for line in f:
                if 'Zero point vibrational energy' in line:
                    zpe = float(line.split()[4]) * 4184  # QChem's ZPE is in kcal/mol, convert to J/mol
                    logging.debug('ZPE is {}'.format(str(zpe)))
        if zpe is not None:
            return zpe
        else:
            raise InputError('Unable to find zero-point energy in QChem output file.')

    def loadScanEnergies(self):
        """
        Extract the optimized energies in J/mol from a QChem log file, e.g. the
        result of a QChem "PES Scan" quantum chemistry calculation.
        """
        Vlist = []
        angle = []
        read = False
        with open(self.path, 'r') as f:
            for line in f:
                if '-----------------' in line:
                    read = False
                if read:
                    values = [float(item) for item in line.split()]
                    angle.append(values[0])
                    Vlist.append(values[1])
                if 'Summary of potential scan:' in line:
                    logging.info('found a sucessfully completed QChem Job')
                    read = True
                elif 'SCF failed to converge' in line:
                    raise InputError('QChem Job did not sucessfully complete: SCF failed to converge')
        logging.info('   Assuming {0} is the output from a QChem PES scan...'.format(os.path.basename(self.path)))

        Vlist = numpy.array(Vlist, numpy.float64)
        # check to see if the scanlog indicates that one of your reacting species may not be the lowest energy conformer
        check_conformer_energy(Vlist, self.path)

        # Adjust energies to be relative to minimum energy conformer
        # Also convert units from Hartree/particle to J/mol
        Vlist -= numpy.min(Vlist)
        Vlist *= constants.E_h * constants.Na
        angle = numpy.arange(0.0, 2 * math.pi + 0.00001, 2 * math.pi / (len(Vlist) - 1), numpy.float64)
        return Vlist, angle

    def loadNegativeFrequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        frequency = 0
        with open(self.path, 'r') as f:
            for line in f:
                # Read imaginary frequency
                if ' Frequency:' in line:
                    frequency = float((line.split()[1]))
                    break
        # Make sure the frequency is imaginary:
        if frequency < 0:
            return frequency
        else:
            raise InputError('Unable to find imaginary frequency in QChem output file {0}'.format(self.path))
