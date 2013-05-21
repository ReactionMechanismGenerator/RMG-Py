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
This module provides the :class:`StatMechJob` class, which represents a
statistical mechanics job used to compute and save the statistical mechanics
information for a single species or transition state.
"""

import os.path
import math
import numpy
import logging

import rmgpy.constants as constants
from rmgpy.cantherm.output import prettify
from rmgpy.cantherm.gaussian import GaussianLog
from rmgpy.species import TransitionState
from rmgpy.statmech import *

################################################################################

class InputError(Exception):
    """
    An exception raised when parsing an input file for a conformer. Pass a
    string describing the error.
    """
    pass

################################################################################

class ScanLog:
    """
    Represent a text file containing a table of angles and corresponding
    scan energies.
    """

    angleFactors = {
        'radians': 1.0,
        'rad': 1.0,
        'degrees': 180.0 / math.pi,
        'deg': 180.0 / math.pi,
    }
    energyFactors = {
        'J/mol': 1.0,
        'kJ/mol': 1.0/1000.,
        'cal/mol': 1.0/4.184,
        'kcal/mol': 1.0/4184.,
        'cm^-1': 1.0/(constants.h * constants.c * 100. * constants.Na),
        'hartree': 1.0/(constants.E_h * constants.Na),
    }
        
    def __init__(self, path):
        self.path = path

    def load(self):
        """
        Load the scan energies from the file. Returns arrays containing the
        angles (in radians) and energies (in J/mol).
        """
        angles = []; energies = []
        angleUnits = None; energyUnits = None
        angleFactor = None; energyFactor = None
        
        with open(self.path, 'r') as stream:
            for line in stream:
                line = line.strip()
                if line == '': continue
                
                tokens = line.split()
                if angleUnits is None or energyUnits is None:
                    angleUnits = tokens[1][1:-1]
                    energyUnits = tokens[3][1:-1]
                    
                    try:
                        angleFactor = ScanLog.angleFactors[angleUnits]
                    except KeyError:
                        raise ValueError('Invalid angle units {0!r}.'.format(angleUnits))
                    try:
                        energyFactor = ScanLog.energyFactors[energyUnits]
                    except KeyError:
                        raise ValueError('Invalid energy units {0!r}.'.format(energyUnits))
            
                else:
                    angles.append(float(tokens[0]) / angleFactor)
                    energies.append(float(tokens[1]) / energyFactor)
        
        angles = numpy.array(angles)
        energies = numpy.array(energies)
        energies -= energies[0]
        
        return angles, energies
        
    def save(self, angles, energies, angleUnits='radians', energyUnits='kJ/mol'):
        """
        Save the scan energies to the file using the given `angles` in radians
        and corresponding energies `energies` in J/mol. The file is created to
        use the given `angleUnits` for angles and `energyUnits` for energies.
        """
        assert len(angles) == len(energies)
        
        try:
            angleFactor = ScanLog.angleFactors[angleUnits]
        except KeyError:
            raise ValueError('Invalid angle units {0!r}.'.format(angleUnits))
        try:
            energyFactor = ScanLog.energyFactors[energyUnits]
        except KeyError:
            raise ValueError('Invalid energy units {0!r}.'.format(energyUnits))
        
        with open(self.path, 'w') as stream:
            stream.write('{0:>24} {1:>24}\n'.format(
                'Angle ({0})'.format(angleUnits),
                'Energy ({0})'.format(energyUnits),
            ))
            for angle, energy in zip(angles, energies):
                stream.write('{0:23.10f} {1:23.10f}\n'.format(angle * angleFactor, energy * energyFactor))

################################################################################

def hinderedRotor(scanLog, pivots, top, symmetry):
    return [scanLog, pivots, top, symmetry]

class StatMechJob:
    """
    A representation of a CanTherm statistical mechanics job. This job is used
    to compute and save the statistical mechanics information for a single
    species or transition state.
    """
    
    def __init__(self, species, path):
        self.species = species
        self.path = path
        self.modelChemistry = ''
        self.frequencyScaleFactor = 1.0
        self.includeHinderedRotors = True
        self.applyBondEnergyCorrections = True
    
    def execute(self, outputFile=None, plot=False):
        """
        Execute the statistical mechanics job, saving the results to the
        given `outputFile` on disk.
        """
        self.load()
        if outputFile is not None:
            self.save(outputFile)
    
    def load(self):
        """
        Load the statistical mechanics parameters for each conformer from
        the associated files on disk. Creates :class:`Conformer` objects for
        each conformer and appends them to the list of confomers on the
        species object.
        """
        logging.info('Loading statistical mechanics parameters for {0}...'.format(self.species.label))
        
        path = self.path
        
        TS = isinstance(self.species, TransitionState)
    
        global_context = {
            '__builtins__': None,
        }
        local_context = {
            '__builtins__': None,
            'True': True,
            'False': False,
            'HinderedRotor': hinderedRotor,
            # File formats
            'GaussianLog': GaussianLog,
            'ScanLog': ScanLog,
        }
    
        directory = os.path.abspath(os.path.dirname(path))
    
        with open(path, 'r') as f:
            try:
                exec f in global_context, local_context
            except (NameError, TypeError, SyntaxError), e:
                logging.error('The species file {0} was invalid:'.format(path))
                raise
        
        try:
            atoms = local_context['atoms']
        except KeyError:
            raise InputError('Required attribute "atoms" not found in species file {0!r}.'.format(path))
        
        try:
            bonds = local_context['bonds']
        except KeyError:
            bonds = {}
            
        try:
            linear = local_context['linear']
        except KeyError:
            raise InputError('Required attribute "linear" not found in species file {0!r}.'.format(path))
        
        try:
            externalSymmetry = local_context['externalSymmetry']
        except KeyError:
            raise InputError('Required attribute "externalSymmetry" not found in species file {0!r}.'.format(path))
        
        try:
            spinMultiplicity = local_context['spinMultiplicity']
        except KeyError:
            raise InputError('Required attribute "spinMultiplicity" not found in species file {0!r}.'.format(path))
       
        try:
            opticalIsomers = local_context['opticalIsomers']
        except KeyError:
            raise InputError('Required attribute "opticalIsomers" not found in species file {0!r}.'.format(path))
        
        try:
            energy = local_context['energy']
        except KeyError:
            raise InputError('Required attribute "energy" not found in species file {0!r}.'.format(path))
        if isinstance(energy, dict):
            try:
                energy = energy[self.modelChemistry]
            except KeyError:
                raise InputError('Model chemistry {0!r} not found in from dictionary of energy values in species file {1!r}.'.format(self.modelChemistry, path))
        if isinstance(energy, GaussianLog):
            energyLog = energy; E0 = None
            energyLog.path = os.path.join(directory, energyLog.path)
        elif isinstance(energy, float):
            energyLog = None; E0 = energy
        
        try:
            geomLog = local_context['geometry']
        except KeyError:
            raise InputError('Required attribute "geometry" not found in species file {0!r}.'.format(path))
        geomLog.path = os.path.join(directory, geomLog.path)
    
        try:
            statmechLog = local_context['frequencies']
        except KeyError:
            raise InputError('Required attribute "frequencies" not found in species file {0!r}.'.format(path))
        statmechLog.path = os.path.join(directory, statmechLog.path)
        
        if 'frequencyScaleFactor' in local_context:
            logging.warning('Ignoring frequency scale factor in species file {0!r}.'.format(path))
        
        try:
            rotors = local_context['rotors']
        except KeyError:
            rotors = []
        
        # But don't consider hindered rotors if flag is not set
        if not self.includeHinderedRotors:
            rotors = []
        
        logging.debug('    Reading molecular degrees of freedom...')
        conformer = statmechLog.loadConformer(symmetry=externalSymmetry, spinMultiplicity=spinMultiplicity, opticalIsomers=opticalIsomers)
        
        logging.debug('    Reading optimized geometry...')
        coordinates, number, mass = geomLog.loadGeometry()
        conformer.coordinates = (coordinates,"angstroms") 
        conformer.number = number
        conformer.mass = (mass,"amu")
        
        logging.debug('    Reading energy...')
        if E0 is None:
            E0 = energyLog.loadEnergy()
        else:
            E0 = E0 * constants.E_h * constants.Na         # Hartree/particle to J/mol
        E0 = applyEnergyCorrections(E0, self.modelChemistry, atoms, bonds if self.applyBondEnergyCorrections else {})
        ZPE = statmechLog.loadZeroPointEnergy() * self.frequencyScaleFactor
        E0 += ZPE
        logging.debug('         ZPE (0 K) = {0:g} kcal/mol'.format(ZPE / 4184.))
        logging.debug('         E0 (0 K) = {0:g} kcal/mol'.format(E0 / 4184.))
       
        conformer.E0 = (E0*0.001,"kJ/mol")
        
        # If loading a transition state, also read the imaginary frequency
        if TS:
            self.species.frequency = (statmechLog.loadNegativeFrequency() * self.frequencyScaleFactor, "cm^-1")

        # Read and fit the 1D hindered rotors if applicable
        # If rotors are found, the vibrational frequencies are also
        # recomputed with the torsional modes removed
        F = statmechLog.loadForceConstantMatrix()
        if F is not None and len(mass) > 1 and len(rotors) > 0:
            
            logging.debug('    Fitting {0} hindered rotors...'.format(len(rotors)))
            rotorCount = 0
            for scanLog, pivots, top, symmetry in rotors:
                
                # Load the hindered rotor scan energies
                if isinstance(scanLog, GaussianLog):
                    scanLog.path = os.path.join(directory, scanLog.path)
                    Vlist, angle = scanLog.loadScanEnergies()
                    scanLogOutput = ScanLog(os.path.join(directory, '{0}_rotor_{1}.txt'.format(self.species.label, rotorCount+1)))
                    scanLogOutput.save(angle, Vlist)
                elif isinstance(scanLog, ScanLog):
                    scanLog.path = os.path.join(directory, scanLog.path)
                    angle, Vlist = scanLog.load()
                else:
                    raise Exception('Invalid log file type {0} for scan log.'.format(scanLog.__class__))
                    
                inertia = conformer.getInternalReducedMomentOfInertia(pivots, top) * constants.Na * 1e23
                
                cosineRotor = HinderedRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                cosineRotor.fitCosinePotentialToData(angle, Vlist)
                fourierRotor = HinderedRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                fourierRotor.fitFourierPotentialToData(angle, Vlist)
                
                Vlist_cosine = numpy.zeros_like(angle)
                Vlist_fourier = numpy.zeros_like(angle)
                for i in range(angle.shape[0]):
                    Vlist_cosine[i] = cosineRotor.getPotential(angle[i])
                    Vlist_fourier[i] = fourierRotor.getPotential(angle[i])
                
                rms_cosine = numpy.sqrt(numpy.sum((Vlist_cosine - Vlist) * (Vlist_cosine - Vlist)) / (len(Vlist) - 1)) / 4184.
                rms_fourier = numpy.sqrt(numpy.sum((Vlist_fourier - Vlist) * (Vlist_fourier - Vlist))/ (len(Vlist) - 1)) / 4184.
                
                # Keep the rotor with the most accurate potential
                rotor = cosineRotor if rms_cosine < rms_fourier else fourierRotor
                # However, keep the cosine rotor if it is accurate enough, the
                # fourier rotor is not significantly more accurate, and the cosine
                # rotor has the correct symmetry 
                if rms_cosine < 0.05 and rms_cosine / rms_fourier < 2.0 and rms_cosine / rms_fourier < 4.0 and symmetry == cosineRotor.symmetry:
                    rotor = cosineRotor
                
                conformer.modes.append(rotor)
                
                self.plotHinderedRotor(angle, Vlist, cosineRotor, fourierRotor, rotor, rotorCount, directory)
                
                rotorCount += 1
                       
            logging.debug('    Determining frequencies from reduced force constant matrix...')
            frequencies = numpy.array(projectRotors(conformer, F, rotors, linear, TS))
            
        elif len(conformer.modes) > 2:
            frequencies = conformer.modes[2].frequencies.value_si
            rotors = numpy.array([])
        else:
            frequencies = numpy.array([])
            rotors = numpy.array([])
    
        for mode in conformer.modes:
            if isinstance(mode, HarmonicOscillator):
                mode.frequencies = (frequencies * self.frequencyScaleFactor,"cm^-1")
                
        self.species.conformer = conformer
    
    def save(self, outputFile):
        """
        Save the results of the statistical mechanics job to the file located
        at `path` on disk.
        """
        
        logging.info('Saving statistical mechanics parameters for {0}...'.format(self.species.label))
        f = open(outputFile, 'a')
    
        numbers = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 14: 'Si', 15: 'P', 16: 'S'}
        
        conformer = self.species.conformer
            
        coordinates = conformer.coordinates.value_si * 1e10
        number = conformer.number.value_si
        
        f.write('# Coordinates for {0} (angstroms):\n'.format(self.species.label))
        for i in range(coordinates.shape[0]):
            x = coordinates[i,0] - coordinates[0,0]
            y = coordinates[i,1] - coordinates[0,1]
            z = coordinates[i,2] - coordinates[0,2]
            f.write('#   {0} {1:9.4f} {2:9.4f} {3:9.4f}\n'.format(numbers[number[i]], x, y, z))
        
        string = 'conformer(label={0!r}, E0={1!r}, modes={2!r}, spinMultiplicity={3:d}, opticalIsomers={4:d}'.format(
            self.species.label, 
            conformer.E0,
            conformer.modes,
            conformer.spinMultiplicity,
            conformer.opticalIsomers,
        )
        try:
            string += ', frequency={0!r}'.format(self.species.frequency)
        except AttributeError: pass
        string += ')'
        
        f.write('{0}\n\n'.format(prettify(string)))
        
        f.close()

    def plotHinderedRotor(self, angle, Vlist, cosineRotor, fourierRotor, rotor, rotorIndex, directory):
        """
        Plot the potential for the rotor, along with its cosine and Fourier
        series potential fits. The plot is saved to a set of files of the form
        ``hindered_rotor_1.pdf``.
        """
        try:
            import pylab
        except ImportError:
            return
        
        phi = numpy.arange(0, 6.3, 0.02, numpy.float64)
        Vlist_cosine = numpy.zeros_like(phi)
        Vlist_fourier = numpy.zeros_like(phi)
        for i in range(phi.shape[0]):
            Vlist_cosine[i] = cosineRotor.getPotential(phi[i])
            Vlist_fourier[i] = fourierRotor.getPotential(phi[i])
        
        fig = pylab.figure(figsize=(6,5))
        pylab.plot(angle, Vlist / 4184., 'ok')
        linespec = '-r' if rotor is cosineRotor else '--r'
        pylab.plot(phi, Vlist_cosine / 4184., linespec)
        linespec = '-b' if rotor is fourierRotor else '--b'
        pylab.plot(phi, Vlist_fourier / 4184., linespec)
        pylab.legend(['scan', 'cosine', 'fourier'], loc=1)
        pylab.xlim(0, 2*constants.pi)
        pylab.xlabel('Angle')
        pylab.ylabel('Potential (kcal/mol)')
        pylab.title('{0} hindered rotor #{1:d}'.format(self.species.label, rotorIndex+1))
        
        axes = fig.get_axes()[0]
        axes.set_xticks([float(j*constants.pi/4) for j in range(0,9)])
        axes.set_xticks([float(j*constants.pi/8) for j in range(0,17)], minor=True)
        axes.set_xticklabels(['$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'])
        
        pylab.savefig(os.path.join(directory, '{0}_rotor_{1:d}.pdf'.format(self.species.label, rotorIndex+1)))
        pylab.close()

################################################################################

def applyEnergyCorrections(E0, modelChemistry, atoms, bonds):
    """
    Given an energy `E0` in J/mol as read from the output of a quantum chemistry
    calculation at a given `modelChemistry`, adjust the energy such that it
    is consistent with the normal gas-phase reference states. `atoms` is a
    dictionary associating element symbols with the number of that element in
    the molecule. `bonds` is a dictionary associating bond types with the number
    of that bond in the molecule.
    """
    
    # Spin orbit correction (SOC) in Hartrees
    # Values taken from note 22 of http://jcp.aip.org/resource/1/jcpsa6/v109/i24/p10570_s1 and converted to hartrees
    # Values in millihartree are also available (with fewer significant figures) from http://jcp.aip.org/resource/1/jcpsa6/v106/i3/p1063_s1
    SOC = {'H':0.0, 'N':0.0, 'O': -0.000355, 'C': -0.000135, 'S':  -0.000893, 'P': 0.0} 
    
    # Step 1: Reference all energies to a model chemistry-independent basis
    # by subtracting out that model chemistry's atomic energies
    # Note: If your model chemistry does not include spin orbit coupling, you should add the corrections to the energies here
    if modelChemistry == 'CBS-QB3':
        atomEnergies = {'H':-0.499818 , 'N':-54.520543, 'O':-74.987624, 'C':-37.785385, 'P':-340.817186, 'S': -397.657360}
    elif modelChemistry == 'G3':
        atomEnergies = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432, 'S': -397.961110}
    elif modelChemistry == 'Klip_1':
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'O':-75.00915718 + SOC['O'], 'C':-37.79249556 + SOC['C']}
    elif modelChemistry == 'Klip_2':
        #Klip QCI(tz,qz)
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'O':-75.00692746 + SOC['O'], 'C':-37.79044863 + SOC['C']}
    elif modelChemistry == 'Klip_2_cc':
        #Klip CCSD(T)(tz,qz)
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'O':-75.00681155 + SOC['O'], 'C':-37.79029443 + SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVTZ-F12':
        # 'CCSD(T)-F12/cc-pVTZ-F12' calculated by CCLass
        atomEnergies = {'H':-0.49994557 + SOC['H'], 'N':-54.43186873 + SOC['N'], 'O':-74.92259120 + SOC['O'], 
                        'C':-37.73766692 + SOC['C'], 'P':-340.77098640 + SOC['P'], 'S': -397.62496049 + SOC['S']}
    else:
        logging.warning('Unknown model chemistry "{0}"; not applying energy corrections.'.format(modelChemistry))
        return E0
    for symbol, count in atoms.items():
        if symbol in atomEnergies: E0 -= count * atomEnergies[symbol] * 4.35974394e-18 * constants.Na
        else:
            logging.warning('Ignored unknown atom type "{0}".'.format(symbol))
    
    # Step 2: Atom energy corrections to reach gas-phase reference state
    # Experimental enthalpy of formation at 0 K 
    # See Gaussian thermo whitepaper at http://www.gaussian.com/g_whitepap/thermo.htm)
    # Note: these values are relatively old and some improvement may be possible by using newer values, particularly for carbon
    # However, care should be taken to ensure that they are compatible with the BAC values (if BACs are used)
    atomHf = {'H': 51.63 , 'N': 112.53 ,'O': 58.99 ,'C': 169.98, 'S': 65.55 }
    # Thermal contribution to enthalpy Hss(298 K) - Hss(0 K) reported by Gaussian thermo whitepaper
    # This will be subtracted from the corresponding value in atomHf to produce an enthalpy used in calculating the enthalpy of formation at 298 K
    atomThermal = {'H': 1.01 , 'N': 1.04, 'O': 1.04 ,'C': 0.25, 'S': 1.05 }
    # Total energy correction used to reach gas-phase reference state
    # Note: Spin orbit coupling no longer included in these energies, since some model chemistries include it automatically
    atomEnergies = {}
    for element in atomHf:
        atomEnergies[element] = atomHf[element] - atomThermal[element]
    for symbol, count in atoms.items():
        if symbol in atomEnergies: E0 += count * atomEnergies[symbol] * 4184.
    
    # Step 3: Bond energy corrections
    bondEnergies = { 'C-H': -0.11, 'C-C': -0.3, 'C=C': -0.08, 'C#C': -0.64,
        'O-H': 0.02, 'C-O': 0.33, 'C=O': 0.55, 'N#N': -2.0, 'O=O': -0.2, 
        'H-H': 1.1, 'C#N': -0.89, 'C-S': 0.43, 'S=O': -0.78 }
    for symbol, count in bonds.items():
        if symbol in bondEnergies: E0 += count * bondEnergies[symbol] * 4184.
        else:
            logging.warning('Ignored unknown bond type {0!r}.'.format(symbol))
    
    return E0

def projectRotors(conformer, F, rotors, linear, TS):
    """
    For a given `conformer` with associated force constant matrix `F`, lists of
    rotor information `rotors`, `pivots`, and `top1`, and the linearity of the
    molecule `linear`, project out the nonvibrational modes from the force
    constant matrix and use this to determine the vibrational frequencies. The
    list of vibrational frequencies is returned in cm^-1.
    """
    
    Nrotors = len(rotors)
    Natoms = len(conformer.mass.value)
    Nvib = 3 * Natoms - (5 if linear else 6) - Nrotors - (1 if (TS) else 0)
    mass = conformer.mass.value_si
    coordinates = conformer.coordinates.value_si
    
    if linear:
        D = numpy.zeros((Natoms*3,5+Nrotors), numpy.float64)
    else:
        D = numpy.zeros((Natoms*3,6+Nrotors), numpy.float64)

    for i in range(Natoms):
        # Projection vectors for translation
        D[3*i+0,0] = 1.0
        D[3*i+1,1] = 1.0
        D[3*i+2,2] = 1.0
        # Projection vectors for [external] rotation
        D[3*i:3*i+3,3] = numpy.array([0, -coordinates[i,2], coordinates[i,1]], numpy.float64)
        D[3*i:3*i+3,4] = numpy.array([coordinates[i,2], 0, -coordinates[i,0]], numpy.float64)
        if not linear:
            D[3*i:3*i+3,5] = numpy.array([-coordinates[i,1], coordinates[i,0], 0], numpy.float64)
    for i, rotor in enumerate(rotors):
        scanLog, pivots, top, symmetry = rotor
        # Determine pivot atom
        if pivots[0] in top: pivot = pivots[0]
        elif pivots[1] in top: pivot = pivots[1]
        else: raise Exception('Could not determine pivot atom.')
        # Projection vectors for internal rotation
        e12 = coordinates[pivots[0]-1,:] - coordinates[pivots[1]-1,:]
        e12 /= numpy.linalg.norm(e12)
        for atom in top:
            e31 = coordinates[atom-1,:] - coordinates[pivot-1,:]
            D[3*(atom-1):3*(atom-1)+3,-Nrotors+i] = numpy.cross(e31, e12)

    # Make sure projection matrix is orthonormal
    import scipy.linalg
    D = scipy.linalg.orth(D)

    # Project out the non-vibrational modes from the force constant matrix
    P = numpy.dot(D, D.transpose())
    I = numpy.identity(Natoms*3, numpy.float64)
    F = numpy.dot(I - P, numpy.dot(F, I - P))

    # Generate mass-weighted force constant matrix
    # This converts the axes to mass-weighted Cartesian axes
    # Units of Fm are J/m^2*kg = 1/s^2
    Fm = F.copy()
    for i in range(Natoms):
        for j in range(Natoms):
            for u in range(3):
                for v in range(3):
                    Fm[3*i+u,3*j+v] /= math.sqrt(mass[i] * mass[j])

    # Get eigenvalues of mass-weighted force constant matrix
    eig, V = numpy.linalg.eigh(Fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation
    return numpy.sqrt(eig[-Nvib:]) / (2 * math.pi * constants.c * 100)
