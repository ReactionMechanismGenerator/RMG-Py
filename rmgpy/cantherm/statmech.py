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
from rmgpy.cantherm.molepro import MoleProLog 
from rmgpy.cantherm.qchem import QchemLog 

from rmgpy.species import TransitionState

from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor
from rmgpy.statmech.conformer import Conformer

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

def hinderedRotor(scanLog, pivots, top, symmetry, fit='best'):
    return [scanLog, pivots, top, symmetry, fit]

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
            'QchemLog': QchemLog,
            'MoleProLog': MoleProLog,
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
        elif isinstance(energy, QchemLog):
            energyLog = energy; E0 = None
            energyLog.path = os.path.join(directory, energyLog.path)
        elif isinstance(energy, MoleProLog):
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
        # The E0 that is read from the log file is without the ZPE and corresponds to E_elec
        if E0 is None:
            E0 = energyLog.loadEnergy(self.frequencyScaleFactor)
        else:
            E0 = E0 * constants.E_h * constants.Na         # Hartree/particle to J/mol
        E0 = applyEnergyCorrections(E0, self.modelChemistry, atoms, bonds if self.applyBondEnergyCorrections else {})
        ZPE = statmechLog.loadZeroPointEnergy() * self.frequencyScaleFactor
        
        # The E0_withZPE at this stage contains the ZPE
        E0_withZPE = E0 + ZPE
        
        logging.debug('         Scaling factor used = {0:g}'.format(self.frequencyScaleFactor))
        logging.debug('         ZPE (0 K) = {0:g} kcal/mol'.format(ZPE / 4184.))
        logging.debug('         E0 (0 K) = {0:g} kcal/mol'.format(E0_withZPE / 4184.))
       
        conformer.E0 = (E0_withZPE*0.001,"kJ/mol")
        
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
            for scanLog, pivots, top, symmetry, fit in rotors:
                
                # Load the hindered rotor scan energies
                if isinstance(scanLog, GaussianLog):
                    scanLog.path = os.path.join(directory, scanLog.path)
                    Vlist, angle = scanLog.loadScanEnergies()
                    scanLogOutput = ScanLog(os.path.join(directory, '{0}_rotor_{1}.txt'.format(self.species.label, rotorCount+1)))
                    scanLogOutput.save(angle, Vlist)
                elif isinstance(scanLog, QchemLog):
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
                
                if fit=='cosine':
                    rotor=cosineRotor
                elif fit =='fourier':
                    rotor=fourierRotor
                elif fit =='best':
                
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
            
            # The frequencies have changed after projection, hence we need to recompute the ZPE
            # We might need to multiply the scaling factor to the frequencies 
            ZPE = self.getZPEfromfrequencies(frequencies)
            E0_withZPE = E0 + ZPE
            # Reset the E0 of the conformer
            conformer.E0 = (E0_withZPE*0.001,"kJ/mol")

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
        
    def getZPEfromfrequencies(self, frequencies):
                
        ZPE = 0.0
        
        for freq in frequencies:
            if freq > 0.0:
                ZPE += 0.5 * constants.h * freq * 100.0 * constants.c * constants.Na
                
        return ZPE
        
    
    def save(self, outputFile):
        """
        Save the results of the statistical mechanics job to the file located
        at `path` on disk.
        """
        
        logging.info('Saving statistical mechanics parameters for {0}...'.format(self.species.label))
        f = open(outputFile, 'a')
    
        numbers = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl'}
        
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
        atomEnergies = {'H':-0.499818 + SOC['H'], 'N':-54.520543 + SOC['N'], 'O':-74.987624+ SOC['O'], 'C':-37.785385+ SOC['C'], 'P':-340.817186+ SOC['P'], 'S': -397.657360+ SOC['S']}
    elif modelChemistry == 'G3':
        atomEnergies = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432, 'S': -397.961110}
    elif modelChemistry == 'M08SO/MG3S*': # * indicates that the grid size used in the [QChem} electronic 
        #structure calculation utilized 75 radial points and 434 angular points 
        #(i.e,, this is specified in the $rem section of the [qchem] input file as: XC_GRID 000075000434)
        atomEnergies = {'H':-0.5017321350 + SOC['H'], 'N':-54.5574039365 + SOC['N'], 'O':-75.0382931348+ SOC['O'], 'C':-37.8245648740+ SOC['C'], 'P':-341.2444299005+ SOC['P'], 'S':-398.0940312227+ SOC['S'] }
    elif modelChemistry == 'Klip_1':
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'N':-54.53383153 + SOC['N'], 'O':-75.00935474+ SOC['O'], 'C':-37.79266591+ SOC['C']}
    elif modelChemistry == 'Klip_2':
        #Klip QCI(tz,qz)
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'N':-54.53169400 + SOC['N'], 'O':-75.00714902+ SOC['O'], 'C':-37.79060419+ SOC['C']}
    elif modelChemistry == 'Klip_3':
        #Klip QCI(dz,tz)
        atomEnergies = {'H':-0.50005578 + SOC['H'], 'N':-54.53128140 + SOC['N'], 'O':-75.00356581+ SOC['O'], 'C':-37.79025175+ SOC['C']}

    elif modelChemistry == 'Klip_2_cc':
        #Klip CCSD(T)(tz,qz)
        atomEnergies = {'H':-0.50003976 + SOC['H'], 'O':-75.00681155+ SOC['O'], 'C':-37.79029443+ SOC['C']}

    elif modelChemistry == 'CCSD(T)-F12/cc-pVDZ-F12_H-TZ':
        atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.526406291655 + SOC['N'], 'O':-74.995458316117+ SOC['O'], 'C':-37.788203485235+ SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVDZ-F12_H-QZ':
        atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.526406291655 + SOC['N'], 'O':-74.995458316117+ SOC['O'], 'C':-37.788203485235+ SOC['C']}
    
    # We are assuming that SOC is included in the Bond Energy Corrections  
    elif modelChemistry == 'CCSD(T)-F12/cc-pVDZ-F12':
#        atomEnergies = {'H':-0.499811124128, 'N':-54.526406291655, 'O':-74.995458316117, 'C':-37.788203485235}
        atomEnergies = {'H':-0.499811124128, 'N':-54.526406291655, 'O':-74.995458316117, 'C':-37.788203485235, 'S':-397.663040369707}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVTZ-F12':
        atomEnergies = {'H':-0.499946213243, 'N':-54.53000909621, 'O':-75.004127673424, 'C':-37.789862146471, 'S':-397.675447487865}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVQZ-F12':
        atomEnergies = {'H':-0.499994558325, 'N':-54.530515226371, 'O':-75.005600062003, 'C':-37.789961656228, 'S':-397.676719774973}
    elif modelChemistry == 'CCSD(T)-F12/cc-pCVDZ-F12':
        atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.582137180344 + SOC['N'], 'O':-75.053045547421 + SOC['O'], 'C':-37.840869118707+ SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/cc-pCVTZ-F12':
        atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.588545831900 + SOC['N'], 'O':-75.065995072347 + SOC['O'], 'C':-37.844662139972+ SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/cc-pCVQZ-F12':
        atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.589137594139+ SOC['N'], 'O':-75.067412234737+ SOC['O'], 'C':-37.844893820561+ SOC['C']}

    elif modelChemistry == 'CCSD(T)-F12/aug-cc-pVDZ':
        atomEnergies = {'H':-0.499459066131 + SOC['H'], 'N':-54.524279516472 + SOC['N'], 'O':-74.992097308083+ SOC['O'], 'C':-37.786694171716+ SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/aug-cc-pVTZ':
        atomEnergies = {'H':-0.499844820798 + SOC['H'], 'N':-54.527419359906 + SOC['N'], 'O':-75.000001429806+ SOC['O'], 'C':-37.788504810868+ SOC['C']}
    elif modelChemistry == 'CCSD(T)-F12/aug-cc-pVQZ':
        atomEnergies = {'H':-0.499949526073 + SOC['H'], 'N':-54.529569719016 + SOC['N'], 'O':-75.004026586610+ SOC['O'], 'C':-37.789387892348+ SOC['C']}


    elif modelChemistry == 'B-CCSD(T)-F12/cc-pVDZ-F12':
        atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.523269942190 + SOC['N'], 'O':-74.990725918500 + SOC['O'], 'C':-37.785409916465 + SOC['C'], 'S': -397.658155086033 + SOC['S']}
    elif modelChemistry == 'B-CCSD(T)-F12/cc-pVTZ-F12':
        atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.528135889213 + SOC['N'], 'O':-75.001094055506 + SOC['O'], 'C':-37.788233578503 + SOC['C'], 'S':-397.671745425929 + SOC['S']}
    elif modelChemistry == 'B-CCSD(T)-F12/cc-pVQZ-F12':
        atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.529425753163 + SOC['N'], 'O':-75.003820485005 + SOC['O'], 'C':-37.789006506290 + SOC['C'], 'S':-397.674145126931 + SOC['S']}
    elif modelChemistry == 'B-CCSD(T)-F12/cc-pCVDZ-F12':
        atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.578602780288 + SOC['N'], 'O':-75.048064317367+ SOC['O'], 'C':-37.837592033417+ SOC['C']}
    elif modelChemistry == 'B-CCSD(T)-F12/cc-pCVTZ-F12':
        atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.586402551258 + SOC['N'], 'O':-75.062767632757+ SOC['O'], 'C':-37.842729156944+ SOC['C']}
    elif modelChemistry == 'B-CCSD(T)-F12/cc-pCVQZ-F12':
        atomEnergies = {'H':-0.49999456 + SOC['H'], 'N':-54.587781507581 + SOC['N'], 'O':-75.065397706471+ SOC['O'], 'C':-37.843634971592+ SOC['C']}

    elif modelChemistry == 'B-CCSD(T)-F12/aug-cc-pVDZ':
        atomEnergies = {'H':-0.499459066131 + SOC['H'], 'N':-54.520475581942 + SOC['N'], 'O':-74.986992215049+ SOC['O'], 'C':-37.783294495799+ SOC['C']}
    elif modelChemistry == 'B-CCSD(T)-F12/aug-cc-pVTZ':
        atomEnergies = {'H':-0.499844820798 + SOC['H'], 'N':-54.524927371700 + SOC['N'], 'O':-74.996328829705+ SOC['O'], 'C':-37.786320700792+ SOC['C']}
    elif modelChemistry == 'B-CCSD(T)-F12/aug-cc-pVQZ':
        atomEnergies = {'H':-0.499949526073 + SOC['H'], 'N':-54.528189769291 + SOC['N'], 'O':-75.001879610563+ SOC['O'], 'C':-37.788165047059+ SOC['C']}

    elif modelChemistry == 'DFT_G03_b3lyp':
        atomEnergies = {'H':-0.502256981529 + SOC['H'], 'N':-54.6007233648 + SOC['N'], 'O':-75.0898777574+ SOC['O'], 'C':-37.8572666349+ SOC['C']}
    elif modelChemistry == 'DFT_ks_b3lyp':
        atomEnergies = {'H':-0.49785866 + SOC['H'], 'N':-54.45608798 + SOC['N'], 'O':-74.93566254+ SOC['O'], 'C':-37.76119132+ SOC['C']}
    elif modelChemistry == 'DFT_uks_b3lyp':
        atomEnergies = {'H':-0.49785866 + SOC['H'], 'N':-54.45729113 + SOC['N'], 'O':-74.93566254+ SOC['O'], 'C':-37.76119132+ SOC['C']}

    elif modelChemistry == 'MP2_rmp2_pVDZ':
        atomEnergies = {'H':-0.49927840 + SOC['H'], 'N':-54.46141996 + SOC['N'], 'O':-74.89408254+ SOC['O'], 'C':-37.73792713+ SOC['C']}
    elif modelChemistry == 'MP2_rmp2_pVTZ':
        atomEnergies = {'H':-0.49980981 + SOC['H'], 'N':-54.49615972 + SOC['N'], 'O':-74.95506980+ SOC['O'], 'C':-37.75833104+ SOC['C']}
    elif modelChemistry == 'MP2_rmp2_pVQZ':
        atomEnergies = {'H':-0.49994557 + SOC['H'], 'N':-54.50715868 + SOC['N'], 'O':-74.97515364+ SOC['O'], 'C':-37.76533215+ SOC['C']}

    elif modelChemistry == 'CCSD-F12/cc-pVDZ-F12':
        atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.524325513811 + SOC['N'], 'O':-74.992326577897+ SOC['O'], 'C':-37.786213495943+ SOC['C']}

    elif modelChemistry == 'CCSD(T)-F12/cc-pVDZ-F12_noscale':
        atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.526026290887 + SOC['N'], 'O':-74.994751897699+ SOC['O'], 'C':-37.787881871511+ SOC['C']}

    elif modelChemistry == 'G03_PBEPBE_6-311++g_d_p':
        atomEnergies = {'H':-0.499812273282 + SOC['H'], 'N':-54.5289567564 + SOC['N'], 'O':-75.0033596764+ SOC['O'], 'C':-37.7937388736+ SOC['C']}

    elif modelChemistry == 'FCI/cc-pVDZ':
#        atomEnergies = {'C':-37.760717371923}
        atomEnergies = {'C':-37.789527+ SOC['C']}
    elif modelChemistry == 'FCI/cc-pVTZ':
        atomEnergies = {'C':-37.781266669684+ SOC['C']}
    elif modelChemistry == 'FCI/cc-pVQZ':
        atomEnergies = {'C':-37.787052110598+ SOC['C']}
        
    elif modelChemistry in ['BMK/cbsb7', 'BMK/6-311G(2d,d,p)']:
        atomEnergies = {'H':-0.498618853119+ SOC['H'], 'N':-54.5697851544+ SOC['N'], 'O':-75.0515210278+ SOC['O'], 'C':-37.8287310027+ SOC['C'], 'P':-341.167615941+ SOC['P'], 'S': -398.001619915+ SOC['S']}
        
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
    atomHf = {'H': 51.63 , 'N': 112.53 ,'O': 58.99 ,'C': 169.98, 'S': 65.66 }
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
    bondEnergies = {}
    if modelChemistry == 'CCSD(T)-F12/cc-pVDZ-F12':
        bondEnergies = { 'C-H': -0.46, 'C-C': -0.68, 'C=C': -1.90, 'C#C': -3.13,
            'O-H': -0.51, 'C-O': -0.23, 'C=O': -0.69, 'O-O': -0.02, 'N-C': -0.67,
            'N=C': -1.46, 'N#C': -2.79, 'N-O': 0.74, 'N_O': -0.23, 'N=O': -0.51,
            'N-H': -0.69, 'N-N': -0.47, 'N=N': -1.54, 'N#N': -2.05,}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVTZ-F12':
        bondEnergies = { 'C-H': -0.09, 'C-C': -0.27, 'C=C': -1.03, 'C#C': -1.79,
            'O-H': -0.06, 'C-O': 0.14, 'C=O': -0.19, 'O-O': 0.16, 'N-C': -0.18,
            'N=C': -0.41, 'N#C': -1.41, 'N-O': 0.87, 'N_O': -0.09, 'N=O': -0.23,
            'N-H': -0.01, 'N-N': -0.21, 'N=N': -0.44, 'N#N': -0.76,}
    elif modelChemistry == 'CCSD(T)-F12/cc-pVQZ-F12':
        bondEnergies = { 'C-H': -0.08, 'C-C': -0.26, 'C=C': -1.01, 'C#C': -1.66,
            'O-H':  0.07, 'C-O': 0.25, 'C=O': -0.03, 'O-O': 0.26, 'N-C': -0.20,
            'N=C': -0.30, 'N#C': -1.33, 'N-O': 1.01, 'N_O': -0.03, 'N=O': -0.26,
            'N-H':  0.06, 'N-N': -0.23, 'N=N': -0.37, 'N#N': -0.64,}
    
    # BAC corrections from Table IX in http://jcp.aip.org/resource/1/jcpsa6/v109/i24/p10570_s1 for CBS-Q method
    # H-Cl correction from CBS-QB3 enthalpy difference with Gurvich 1989, HF298=-92.31 kJ
    elif modelChemistry == 'CBS-QB3':
        bondEnergies = { 'C-H': -0.11, 'C-C': -0.3, 'C=C': -0.08, 'C#C': -0.64,
            'O-H': 0.02, 'C-O': 0.33, 'C=O': 0.55, 'N#N': -2.0, 'O=O': -0.2, 
            'H-H': 1.1, 'C#N': -0.89, 'C-S': 0.43, 'S=O': -0.78, 'S-H': 0.0, }
    elif modelChemistry in ['B3LYP/cbsb7', 'B3LYP/6-311G(2d,d,p)', 'DFT_G03_b3lyp','B3LYP/6-311+G(3df,2p)']:
        bondEnergies = { 'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
            'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44, 
            'C#N': 0.22, 'C-S': -2.35, 'S=O': -5.19, 'S-H': -0.52, }    
    else:
        
        logging.warning('No bond energy correction found for model chemistry: {0}'.format(modelChemistry))

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
    coordinates = conformer.coordinates.getValue() 


    # Put origin in center of mass
    xm=0.0
    ym=0.0
    zm=0.0
    totmass=0.0
    for i in range(Natoms):
        xm+=mass[i]*coordinates[i,0]
        ym+=mass[i]*coordinates[i,1]
        zm+=mass[i]*coordinates[i,2]
        totmass+=mass[i]
          
    xm/=totmass
    ym/=totmass
    zm/=totmass

    for i in range(Natoms):
        coordinates[i,0]-=xm
        coordinates[i,1]-=ym 
        coordinates[i,2]-=zm
    # Make vector with the root of the mass in amu for each atom
    amass=numpy.sqrt(mass/constants.amu)

    # Rotation matrix
    I=conformer.getMomentOfInertiaTensor()
    PMoI, Ixyz = numpy.linalg.eigh(I)
 
    external=6
    if linear:
        external=5
    
    D = numpy.zeros((Natoms*3,external), numpy.float64)

    P = numpy.zeros((Natoms,3), numpy.float64)

    # Transform the coordinates to the principal axes
    P = numpy.dot(coordinates,Ixyz)

    for i in range(Natoms):
        # Projection vectors for translation
        D[3*i+0,0] = amass[i]
        D[3*i+1,1] = amass[i]
        D[3*i+2,2] = amass[i]

    # Construction of the projection vectors for external rotation
    for i in range(Natoms):
        D[3*i,3] = (P[i,1]*Ixyz[0,2]-P[i,2]*Ixyz[0,1])*amass[i]
        D[3*i+1,3] = (P[i,1]*Ixyz[1,2]-P[i,2]*Ixyz[1,1])*amass[i]
        D[3*i+2,3] = (P[i,1]*Ixyz[2,2]-P[i,2]*Ixyz[2,1])*amass[i]
        D[3*i,4] = (P[i,2]*Ixyz[0,0]-P[i,0]*Ixyz[0,2])*amass[i]
        D[3*i+1,4] = (P[i,2]*Ixyz[1,0]-P[i,0]*Ixyz[1,2])*amass[i]
        D[3*i+2,4] = (P[i,2]*Ixyz[2,0]-P[i,0]*Ixyz[2,2])*amass[i]
        if not linear:
            D[3*i,5] = (P[i,0]*Ixyz[0,1]-P[i,1]*Ixyz[0,0])*amass[i]
            D[3*i+1,5] = (P[i,0]*Ixyz[1,1]-P[i,1]*Ixyz[1,0])*amass[i]
            D[3*i+2,5] = (P[i,0]*Ixyz[2,1]-P[i,1]*Ixyz[2,0])*amass[i]

    # Make sure projection matrix is orthonormal
    import scipy.linalg

    I = numpy.identity(Natoms*3, numpy.float64)

    P = numpy.zeros((Natoms*3,3*Natoms+external), numpy.float64)

    P[:,0:external] = D[:,0:external]
    P[:,external:external+3*Natoms] = I[:,0:3*Natoms]

    for i in range(3*Natoms+external):
        norm=0.0
        for j in range(3*Natoms):
            norm+=P[j,i]*P[j,i]
        for j in range(3*Natoms):
            if (norm>1E-15):
                P[j,i]/=numpy.sqrt(norm)
            else:
                P[j,i]=0.0
        for j in range(i+1,3*Natoms+external):
            proj=0.0
            for k in range(3*Natoms):
                proj+=P[k,i]*P[k,j]
            for k in range(3*Natoms):
                P[k,j]-=proj*P[k,i]

    # Order D, there will be vectors that are 0.0  
    i=0
    while i < 3*Natoms:
        norm=0.0
        for j in range(3*Natoms):
            norm+=P[j,i]*P[j,i]
        if (norm<0.5):
            P[:,i:3*Natoms+external-1] = P[:,i+1:3*Natoms+external]
        else:
            i+=1

    # T is the transformation vector from cartesian to internal coordinates
    T = numpy.zeros((Natoms*3,3*Natoms-external), numpy.float64)

    T[:,0:3*Natoms-external] = P[:,external:3*Natoms]

    # Generate mass-weighted force constant matrix
    # This converts the axes to mass-weighted Cartesian axes
    # Units of Fm are J/m^2*kg = 1/s^2
    Fm = F.copy()
    for i in range(Natoms):
        for j in range(Natoms):
            for u in range(3):
                for v in range(3):
                    Fm[3*i+u,3*j+v] /= math.sqrt(mass[i] * mass[j])

    Fint = numpy.dot(T.T, numpy.dot(Fm,T))

    # Get eigenvalues of internal force constant matrix, V = 3N-6 * 3N-6
    eig, V = numpy.linalg.eigh(Fint)

    logging.debug('Frequencies from internal Hessian')  
    for i in range(3*Natoms-external):
        logging.debug(numpy.sqrt(eig[i])/(2 * math.pi * constants.c * 100))

    # Now we can start thinking about projecting out the internal rotations
    Dint=numpy.zeros((3*Natoms,Nrotors), numpy.float64)

    counter=0
    for i, rotor in enumerate(rotors):
        scanLog, pivots, top, symmetry, fit = rotor
        # Determine pivot atom
        if pivots[0] in top:
            pivot1 = pivots[0]
            pivot2 = pivots[1]
        elif pivots[1] in top:
            pivot1 = pivots[1]
            pivot2 = pivots[0]
        else: raise Exception('Could not determine pivot atom.')
        # Projection vectors for internal rotation
        e12 = coordinates[pivot1-1,:] - coordinates[pivot2-1,:]
        for j in range(Natoms):
            atom=j+1
            if atom in top:
                e31 = coordinates[atom-1,:] - coordinates[pivot1-1,:]
                Dint[3*(atom-1):3*(atom-1)+3,counter] = numpy.cross(e31, e12)*amass[atom-1]
            else:
                e31 = coordinates[atom-1,:] - coordinates[pivot2-1,:]
                Dint[3*(atom-1):3*(atom-1)+3,counter] = numpy.cross(e31, -e12)*amass[atom-1]
        counter+=1

    # Normal modes in mass weighted cartesian coordinates
    Vmw = numpy.dot(T,V)
    eigM = numpy.zeros((3*Natoms-external,3*Natoms-external), numpy.float64)

    for i in range(3*Natoms-external):
        eigM[i,i]=eig[i]
 
    Fm=numpy.dot(Vmw,numpy.dot(eigM,Vmw.T))

    # Internal rotations are not normal modes => project them on the normal modes and orthogonalize
    # Dintproj =  (3N-6) x (3N) x (3N) x (Nrotors)
    Dintproj=numpy.dot(Vmw.T,Dint)    

    # Reconstruct Dint
    for i in range(Nrotors):
        for j in range (3*Natoms):
            Dint[j,i]=0
            for k in range(3*Natoms-external):
                Dint[j,i]+=Dintproj[k,i]*Vmw[j,k]

    # Ortho normalize
    for i in range(Nrotors):
        norm=0.0
        for j in range(3*Natoms):
            norm+=Dint[j,i]*Dint[j,i]
        for j in range(3*Natoms):
            Dint[j,i]/=numpy.sqrt(norm)
        for j in range(i+1,Nrotors):
            proj=0.0
            for k in range (3*Natoms):
                proj+=Dint[k,i]*Dint[k,j]
            for k in range(3*Natoms):
                Dint[k,j]-=proj*Dint[k,i]

    Dintproj=numpy.dot(Vmw.T,Dint)
    Proj = numpy.dot(Dint, Dint.T)
    I = numpy.identity(Natoms*3, numpy.float64)
    Proj = I - Proj 
    Fm=numpy.dot(Proj, numpy.dot(Fm,Proj))
    # Get eigenvalues of mass-weighted force constant matrix
    eig, V = numpy.linalg.eigh(Fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation

    logging.debug('Frequencies from projected Hessian')
    for i in range(3*Natoms):
        logging.debug(numpy.sqrt(eig[i])/(2 * math.pi * constants.c * 100))
        
    return numpy.sqrt(eig[-Nvib:]) / (2 * math.pi * constants.c * 100)

