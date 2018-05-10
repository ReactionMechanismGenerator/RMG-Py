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
This module provides the :class:`StatMechJob` class, which represents a
statistical mechanics job used to compute and save the statistical mechanics
information for a single species or transition state.
"""

import os.path
import math
import numpy
import logging

from rdkit.Chem import GetPeriodicTable

import rmgpy.constants as constants

from rmgpy.cantherm.output import prettify
from rmgpy.cantherm.gaussian import GaussianLog
from rmgpy.cantherm.molpro import MolproLog
from rmgpy.cantherm.qchem import QchemLog 

from rmgpy.species import TransitionState, Species

from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor, FreeRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.exceptions import InputError

# These are the atoms we currently have enthalpies of formation for
atom_num_dict = {1: 'H',
                 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 53: 'I'}

# Use the RDKit periodic table so we can write symbols for not implemented elements
_rdkit_periodic_table = GetPeriodicTable()

################################################################################

class ScanLog(object):
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


def freeRotor(pivots,top,symmetry):
    return [pivots,top,symmetry]


class StatMechJob(object):
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
        self.applyAtomEnergyCorrections = True
        self.applyBondEnergyCorrections = True
        self.atomEnergies = None
    
    def execute(self, outputFile=None, plot=False):
        """
        Execute the statistical mechanics job, saving the results to the
        given `outputFile` on disk.
        """
        self.load()
        if outputFile is not None:
            self.save(outputFile)
        logging.debug('Finished statmech job for species {0}.'.format(self.species))
        logging.debug(repr(self.species))
    
    def load(self):
        """
        Load the statistical mechanics parameters for each conformer from
        the associated files on disk. Creates :class:`Conformer` objects for
        each conformer and appends them to the list of conformers on the
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
            'FreeRotor': freeRotor,
            # File formats
            'GaussianLog': GaussianLog,
            'QchemLog': QchemLog,
            'MolproLog': MolproLog,
            'ScanLog': ScanLog,
            'Log': Log
        }
    
        directory = os.path.abspath(os.path.dirname(path))
    
        with open(path, 'r') as f:
            try:
                exec f in global_context, local_context
            except (NameError, TypeError, SyntaxError), e:
                logging.error('The species file {0} was invalid:'.format(path))
                raise
        
        try:
            bonds = local_context['bonds']
        except KeyError:
            bonds = {}
            
        try:
            linear = local_context['linear']
            symfromlog = False
        except KeyError:
            externalSymmetry = None
            symfromlog = True

        try:
            externalSymmetry = local_context['externalSymmetry']
            symfromlog = False
        except KeyError:
            externalSymmetry = None
            symfromlog = True
            
        try:
            spinMultiplicity = local_context['spinMultiplicity']
        except KeyError:
            spinMultiplicity = 0
       
        try:
            opticalIsomers = local_context['opticalIsomers']
        except KeyError:
            raise InputError('Required attribute "opticalIsomers" not found in species file {0!r}.'.format(path))
        
        try:
            energy = local_context['energy']
        except KeyError:
            raise InputError('Required attribute "energy" not found in species file {0!r}.'.format(path))
        if isinstance(energy, dict):
            energy = {k.lower(): v for k, v in energy.items()}  # Make model chemistries lower-case
            try:
                energy = energy[self.modelChemistry]
            except KeyError:
                raise InputError('Model chemistry {0!r} not found in from dictionary of energy values in species file '
                                 '{1!r}.'.format(self.modelChemistry, path))
        if isinstance(energy, Log):
            energy.determine_qm_software(os.path.join(directory, energy.path))
            energy = energy.software_log
        if isinstance(energy, (GaussianLog,QchemLog,MolproLog)):
            energyLog = energy; E0 = None
            energyLog.path = os.path.join(directory, energyLog.path)
        elif isinstance(energy, float):
            energyLog = None; E0 = energy
        
        try:
            geomLog = local_context['geometry']
        except KeyError:
            raise InputError('Required attribute "geometry" not found in species file {0!r}.'.format(path))
        if isinstance(geomLog, Log):
            geomLog.determine_qm_software(os.path.join(directory, geomLog.path))
            geomLog = geomLog.software_log
        else:
            geomLog.path = os.path.join(directory, geomLog.path)
    
        try:
            statmechLog = local_context['frequencies']
        except KeyError:
            raise InputError('Required attribute "frequencies" not found in species file {0!r}.'.format(path))
        if isinstance(statmechLog, Log):
            statmechLog.determine_qm_software(os.path.join(directory, statmechLog.path))
            statmechLog = statmechLog.software_log
        else:
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

        #If hindered/free rotors are included in Statmech job, ensure that the same (freq) log file is used for
        # both the species's optimized geometry and Hessian. This approach guarantees that the geometry and Hessian
        #will be defined in the same Cartesian coordinate system ("Input Orientation", as opposed to "Standard Orientation",
        #or something else). Otherwise, if the geometry and Hessian are read from different log files, it is very easy
        #for them to be defined in different coordinate systems, unless the user is very careful. The current implementation
        #only performs this check for Gaussian logs. If QChem logs are used, only a warning is output reminding the user
        #to ensure the geometry and Hessian are defined in consistent coordinates.
        if len(rotors) > 0:
            if isinstance(statmechLog, GaussianLog):
                if statmechLog.path != geomLog.path:
                    raise InputError('For {0!r}, the geometry log, {1!r}, and frequency log, {2!r}, are not the same.'
                                     ' In order to ensure the geometry and Hessian of {0!r} are defined in consistent coordinate systems'
                                     ' for hindered/free rotor projection, either use the frequency log for both geometry and frequency,'
                                     ' or remove rotors.'.format(self.species.label,geomLog.path,statmechLog.path))
            elif isinstance(statmechLog, QchemLog):
                    logging.warning('Qchem log will be used for Hessian of {0!r}. '
                                    'Please verify that the geometry and Hessian of {0!r} are defined in the same coordinate system'.format(self.species.label))

        logging.debug('    Reading molecular degrees of freedom...')
        conformer = statmechLog.loadConformer(symmetry=externalSymmetry, spinMultiplicity=spinMultiplicity,
                                              opticalIsomers=opticalIsomers, symfromlog=symfromlog,
                                              label=self.species.label)

        if conformer.spinMultiplicity == 0:
            raise ValueError("Could not read spin multiplicity from log file {0},\n"
                             "please specify the multiplicity in the input file.".format(self.path))

        logging.debug('    Reading optimized geometry...')
        coordinates, number, mass = geomLog.loadGeometry()

        # Infer atoms from geometry
        atoms = {}
        for atom_num in number:
            try:
                symbol = atom_num_dict[atom_num]
            except KeyError:
                raise Exception(
                    'Element {} is not yet supported.'.format(_rdkit_periodic_table.GetElementSymbol(atom_num))
                )
            atoms[symbol] = atoms.get(symbol, 0) + 1

        # Save atoms for use in writing thermo output
        if isinstance(self.species, Species):
            self.species.props['elementCounts'] = atoms

        conformer.coordinates = (coordinates,"angstroms")
        conformer.number = number
        conformer.mass = (mass,"amu")

        logging.debug('    Reading energy...')
        # The E0 that is read from the log file is without the ZPE and corresponds to E_elec
        if E0 is None:
            E0 = energyLog.loadEnergy(self.frequencyScaleFactor)
        else:
            E0 = E0 * constants.E_h * constants.Na         # Hartree/particle to J/mol
        if not self.applyAtomEnergyCorrections:
            logging.warning('Atom corrections are not being used. Do not trust energies and thermo.')
        E0 = applyEnergyCorrections(E0,
                                    self.modelChemistry,
                                    atoms,
                                    bonds,
                                    atomEnergies=self.atomEnergies,
                                    applyAtomEnergyCorrections=self.applyAtomEnergyCorrections,
                                    applyBondEnergyCorrections=self.applyBondEnergyCorrections)
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
            for q in rotors:
                if len(q) == 3:
                    pivots, top, symmetry = q
                    inertia = conformer.getInternalReducedMomentOfInertia(pivots, top) * constants.Na * 1e23
                    rotor = FreeRotor(inertia=(inertia,"amu*angstrom^2"),symmetry=symmetry)
                    conformer.modes.append(rotor)
                    rotorCount += 1
                elif len(q) == 5:
                    scanLog, pivots, top, symmetry, fit  = q
                    # Load the hindered rotor scan energies
                    if isinstance(scanLog, Log):
                        scanLog.determine_qm_software(os.path.join(directory, scanLog.path))
                        scanLog = scanLog.software_log
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
                        rotorCount += 1
                        conformer.modes.append(rotor)
                    elif fit =='fourier':
                        rotor=fourierRotor
                        rotorCount += 1
                        conformer.modes.append(rotor)
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

        elif len(conformer.modes) > 2:
            if len(rotors) > 0:
                logging.warn('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for Gaussian to generate the force constant matrix, if running Molpro include keyword print, hessian')
            frequencies = conformer.modes[2].frequencies.value_si
            rotors = numpy.array([])
        else:
            if len(rotors) > 0:
                logging.warn('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for Gaussian to generate the force constant matrix, if running Molpro include keyword print, hessian')
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

        conformer = self.species.conformer
            
        coordinates = conformer.coordinates.value_si * 1e10
        number = conformer.number.value_si
        
        f.write('# Coordinates for {0} in Input Orientation (angstroms):\n'.format(self.species.label))
        for i in range(coordinates.shape[0]):
            x = coordinates[i,0]
            y = coordinates[i,1]
            z = coordinates[i,2]
            f.write('#   {0} {1:9.4f} {2:9.4f} {3:9.4f}\n'.format(atom_num_dict[number[i]], x, y, z))
        
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

def applyEnergyCorrections(E0, modelChemistry, atoms, bonds,
                           atomEnergies=None, applyAtomEnergyCorrections=True, applyBondEnergyCorrections=False):
    """
    Given an energy `E0` in J/mol as read from the output of a quantum chemistry
    calculation at a given `modelChemistry`, adjust the energy such that it
    is consistent with the normal gas-phase reference states. `atoms` is a
    dictionary associating element symbols with the number of that element in
    the molecule. The atom energies are in Hartrees, which are from single
    atom calculations using corresponding model chemistries. 
    
    The assumption for the multiplicity of each atom is: 
    H singlet, C triplet, O triplet, N quartet, S triplet, P quartet.

    `bonds` is a dictionary associating bond types with the number
    of that bond in the molecule.
    """

    if applyAtomEnergyCorrections:
        # Spin orbit correction (SOC) in Hartrees
        # Values taken from ref 22 of http://dx.doi.org/10.1063/1.477794 and converted to hartrees
        # Values in millihartree are also available (with fewer significant figures) from table VII of http://dx.doi.org/10.1063/1.473182
        # Iodine SOC calculated as a weighted average of the electronic spin splittings of the lowest energy state. The splittings are
        # obtained from Huber, K.P.; Herzberg, G., Molecular Spectra and Molecular Structure. IV. Constants of Diatomic Molecules, Van Nostrand Reinhold Co., 1979
        SOC = {'H':0.0, 'N':0.0, 'O': -0.000355, 'C': -0.000135, 'S':  -0.000893, 'P': 0.0, 'I':-0.011547226,}

        # Step 1: Reference all energies to a model chemistry-independent basis
        # by subtracting out that model chemistry's atomic energies
        # All model chemistries here should be lower-case because the user input is changed to lower-case
        if atomEnergies is None:
            # Note: If your model chemistry does not include spin orbit coupling, you should add the corrections to the energies here
            if modelChemistry == 'cbs-qb3':
                atomEnergies = {'H':-0.499818 + SOC['H'], 'N':-54.520543 + SOC['N'], 'O':-74.987624+ SOC['O'], 'C':-37.785385+ SOC['C'], 'P':-340.817186+ SOC['P'], 'S': -397.657360+ SOC['S']}
            elif modelChemistry == 'm06-2x/cc-pvtz':
                atomEnergies = {'H':-0.498135 + SOC['H'], 'N':-54.586780 + SOC['N'], 'O':-75.064242+ SOC['O'], 'C':-37.842468+ SOC['C'], 'P':-341.246985+ SOC['P'], 'S': -398.101240+ SOC['S']}
            elif modelChemistry == 'g3':
                atomEnergies = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432, 'S': -397.961110}
            elif modelChemistry == 'm08so/mg3s*': # * indicates that the grid size used in the [QChem] electronic
                #structure calculation utilized 75 radial points and 434 angular points
                #(i.e,, this is specified in the $rem section of the [qchem] input file as: XC_GRID 000075000434)
                atomEnergies = {'H':-0.5017321350 + SOC['H'], 'N':-54.5574039365 + SOC['N'], 'O':-75.0382931348+ SOC['O'], 'C':-37.8245648740+ SOC['C'], 'P':-341.2444299005+ SOC['P'], 'S':-398.0940312227+ SOC['S'] }
            elif modelChemistry == 'klip_1':
                atomEnergies = {'H':-0.50003976 + SOC['H'], 'N':-54.53383153 + SOC['N'], 'O':-75.00935474+ SOC['O'], 'C':-37.79266591+ SOC['C']}
            elif modelChemistry == 'klip_2':
                #Klip QCI(tz,qz)
                atomEnergies = {'H':-0.50003976 + SOC['H'], 'N':-54.53169400 + SOC['N'], 'O':-75.00714902+ SOC['O'], 'C':-37.79060419+ SOC['C']}
            elif modelChemistry == 'klip_3':
                #Klip QCI(dz,tz)
                atomEnergies = {'H':-0.50005578 + SOC['H'], 'N':-54.53128140 + SOC['N'], 'O':-75.00356581+ SOC['O'], 'C':-37.79025175+ SOC['C']}

            elif modelChemistry == 'klip_2_cc':
                #Klip CCSD(T)(tz,qz)
                atomEnergies = {'H':-0.50003976 + SOC['H'], 'O':-75.00681155+ SOC['O'], 'C':-37.79029443+ SOC['C']}

            elif modelChemistry == 'ccsd(t)-f12/cc-pvdz-f12_h-tz':
                atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.526406291655 + SOC['N'], 'O':-74.995458316117+ SOC['O'], 'C':-37.788203485235+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/cc-pvdz-f12_h-qz':
                atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.526406291655 + SOC['N'], 'O':-74.995458316117+ SOC['O'], 'C':-37.788203485235+ SOC['C']}

            # We are assuming that SOC is included in the Bond Energy Corrections
            elif modelChemistry == 'ccsd(t)-f12/cc-pvdz-f12':
                atomEnergies = {'H':-0.499811124128, 'N':-54.526406291655, 'O':-74.995458316117, 'C':-37.788203485235, 'S':-397.663040369707}
            elif modelChemistry == 'ccsd(t)-f12/cc-pvtz-f12':
                atomEnergies = {'H':-0.499946213243, 'N':-54.53000909621, 'O':-75.004127673424, 'C':-37.789862146471, 'S':-397.675447487865}
            elif modelChemistry == 'ccsd(t)-f12/cc-pvqz-f12':
                atomEnergies = {'H':-0.499994558325, 'N':-54.530515226371, 'O':-75.005600062003, 'C':-37.789961656228, 'S':-397.676719774973}
            elif modelChemistry == 'ccsd(t)-f12/cc-pcvdz-f12':
                atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.582137180344 + SOC['N'], 'O':-75.053045547421 + SOC['O'], 'C':-37.840869118707+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/cc-pcvtz-f12':
                atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.588545831900 + SOC['N'], 'O':-75.065995072347 + SOC['O'], 'C':-37.844662139972+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/cc-pcvqz-f12':
                atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.589137594139+ SOC['N'], 'O':-75.067412234737+ SOC['O'], 'C':-37.844893820561+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/cc-pvtz-f12(-pp)':
                atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.53000909621 + SOC['N'], 'O':-75.004127673424 + SOC['O'], 'C':-37.789862146471 + SOC['C'], 'S':-397.675447487865 + SOC['S'], 'I':-294.81781766 + SOC['I']}
            #ccsd(t)/aug-cc-pvtz(-pp) atomic energies were fit to a set of 8 small molecules: CH4, CH3OH, H2S, H2O, SO2, HI, I2, CH3I
            elif modelChemistry == 'ccsd(t)/aug-cc-pvtz(-pp)':
                atomEnergies = {'H':-0.499821176024 + SOC['H'], 'O':-74.96738492 + SOC['O'], 'C':-37.77385697 + SOC['C'], 'S':-397.6461604 + SOC['S'], 'I':-294.7958443 + SOC['I']}

            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvdz':
                atomEnergies = {'H':-0.499459066131 + SOC['H'], 'N':-54.524279516472 + SOC['N'], 'O':-74.992097308083+ SOC['O'], 'C':-37.786694171716+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvtz':
                atomEnergies = {'H':-0.499844820798 + SOC['H'], 'N':-54.527419359906 + SOC['N'], 'O':-75.000001429806+ SOC['O'], 'C':-37.788504810868+ SOC['C']}
            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvqz':
                atomEnergies = {'H':-0.499949526073 + SOC['H'], 'N':-54.529569719016 + SOC['N'], 'O':-75.004026586610+ SOC['O'], 'C':-37.789387892348+ SOC['C']}


            elif modelChemistry == 'b-ccsd(t)-f12/cc-pvdz-f12':
                atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.523269942190 + SOC['N'], 'O':-74.990725918500 + SOC['O'], 'C':-37.785409916465 + SOC['C'], 'S': -397.658155086033 + SOC['S']}
            elif modelChemistry == 'b-ccsd(t)-f12/cc-pvtz-f12':
                atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.528135889213 + SOC['N'], 'O':-75.001094055506 + SOC['O'], 'C':-37.788233578503 + SOC['C'], 'S':-397.671745425929 + SOC['S']}
            elif modelChemistry == 'b-ccsd(t)-f12/cc-pvqz-f12':
                atomEnergies = {'H':-0.499994558325 + SOC['H'], 'N':-54.529425753163 + SOC['N'], 'O':-75.003820485005 + SOC['O'], 'C':-37.789006506290 + SOC['C'], 'S':-397.674145126931 + SOC['S']}
            elif modelChemistry == 'b-ccsd(t)-f12/cc-pcvdz-f12':
                atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.578602780288 + SOC['N'], 'O':-75.048064317367+ SOC['O'], 'C':-37.837592033417+ SOC['C']}
            elif modelChemistry == 'b-ccsd(t)-f12/cc-pcvtz-f12':
                atomEnergies = {'H':-0.499946213243 + SOC['H'], 'N':-54.586402551258 + SOC['N'], 'O':-75.062767632757+ SOC['O'], 'C':-37.842729156944+ SOC['C']}
            elif modelChemistry == 'b-ccsd(t)-f12/cc-pcvqz-f12':
                atomEnergies = {'H':-0.49999456 + SOC['H'], 'N':-54.587781507581 + SOC['N'], 'O':-75.065397706471+ SOC['O'], 'C':-37.843634971592+ SOC['C']}

            elif modelChemistry == 'b-ccsd(t)-f12/aug-cc-pvdz':
                atomEnergies = {'H':-0.499459066131 + SOC['H'], 'N':-54.520475581942 + SOC['N'], 'O':-74.986992215049+ SOC['O'], 'C':-37.783294495799+ SOC['C']}
            elif modelChemistry == 'b-ccsd(t)-f12/aug-cc-pvtz':
                atomEnergies = {'H':-0.499844820798 + SOC['H'], 'N':-54.524927371700 + SOC['N'], 'O':-74.996328829705+ SOC['O'], 'C':-37.786320700792+ SOC['C']}
            elif modelChemistry == 'b-ccsd(t)-f12/aug-cc-pvqz':
                atomEnergies = {'H':-0.499949526073 + SOC['H'], 'N':-54.528189769291 + SOC['N'], 'O':-75.001879610563+ SOC['O'], 'C':-37.788165047059+ SOC['C']}

            elif modelChemistry == 'mp2_rmp2_pvdz':
                atomEnergies = {'H':-0.49927840 + SOC['H'], 'N':-54.46141996 + SOC['N'], 'O':-74.89408254+ SOC['O'], 'C':-37.73792713+ SOC['C']}
            elif modelChemistry == 'mp2_rmp2_pvtz':
                atomEnergies = {'H':-0.49980981 + SOC['H'], 'N':-54.49615972 + SOC['N'], 'O':-74.95506980+ SOC['O'], 'C':-37.75833104+ SOC['C']}
            elif modelChemistry == 'mp2_rmp2_pvqz':
                atomEnergies = {'H':-0.49994557 + SOC['H'], 'N':-54.50715868 + SOC['N'], 'O':-74.97515364+ SOC['O'], 'C':-37.76533215+ SOC['C']}

            elif modelChemistry == 'ccsd-f12/cc-pvdz-f12':
                atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.524325513811 + SOC['N'], 'O':-74.992326577897+ SOC['O'], 'C':-37.786213495943+ SOC['C']}

            elif modelChemistry == 'ccsd(t)-f12/cc-pvdz-f12_noscale':
                atomEnergies = {'H':-0.499811124128 + SOC['H'], 'N':-54.526026290887 + SOC['N'], 'O':-74.994751897699+ SOC['O'], 'C':-37.787881871511+ SOC['C']}

            elif modelChemistry == 'g03_pbepbe_6-311++g_d_p':
                atomEnergies = {'H':-0.499812273282 + SOC['H'], 'N':-54.5289567564 + SOC['N'], 'O':-75.0033596764+ SOC['O'], 'C':-37.7937388736+ SOC['C']}

            elif modelChemistry == 'fci/cc-pvdz':
                atomEnergies = {'C':-37.789527+ SOC['C']}
            elif modelChemistry == 'fci/cc-pvtz':
                atomEnergies = {'C':-37.781266669684+ SOC['C']}
            elif modelChemistry == 'fci/cc-pvqz':
                atomEnergies = {'C':-37.787052110598+ SOC['C']}

            elif modelChemistry in ['bmk/cbsb7', 'bmk/6-311g(2d,d,p)']:
                atomEnergies = {'H':-0.498618853119+ SOC['H'], 'N':-54.5697851544+ SOC['N'], 'O':-75.0515210278+ SOC['O'], 'C':-37.8287310027+ SOC['C'], 'P':-341.167615941+ SOC['P'], 'S': -398.001619915+ SOC['S']}
            elif modelChemistry == 'b3lyp/6-31g**':  # Fitted to small molecules
                atomEnergies = {'H':-0.500426155, 'C':-37.850331697831, 'O':-75.0535872748806, 'S':-398.100820107242}
            elif modelChemistry == 'b3lyp/6-311+g(3df,2p)':  # Calculated atomic energies
                atomEnergies = {'H':-0.502155915123 + SOC['H'], 'C':-37.8574709934 + SOC['C'], 'N':-54.6007233609 + SOC['N'], 'O':-75.0909131284 + SOC['O'], 'P':-341.281730319 + SOC['P'], 'S':-398.134489850 + SOC['S']}

            else:
                raise Exception('Unknown model chemistry "{}".'.format(modelChemistry))

        for symbol, count in atoms.items():
            if symbol in atomEnergies:
                E0 -= count * atomEnergies[symbol] * 4.35974394e-18 * constants.Na
            else:
                raise Exception(
                    'Unknown element "{}". Turn off atom corrections if only running a kinetics jobs '
                    'or supply a dictionary of atom energies.'.format(symbol)
                )

        # Step 2: Atom energy corrections to reach gas-phase reference state
        # Experimental enthalpy of formation at 0 K
        # See Gaussian thermo whitepaper at http://www.gaussian.com/g_whitepap/thermo.htm)
        # Note: these values are relatively old and some improvement may be possible by using newer values, particularly for carbon
        # However, care should be taken to ensure that they are compatible with the BAC values (if BACs are used)
        # The enthalpies listed here should correspond to the allowed elements in atom_num_dict
        # Iodine value is from Cox, J. D., Wagman, D. D., and Medvedev, V. A., CODATA Key Values for Thermodynamics, Hemisphere Publishing Corp., New York, 1989.
        atomHf = {'H': 51.63,
                  'Li': 37.69, 'Be': 76.48, 'B': 136.2, 'C': 169.98, 'N': 112.53, 'O': 58.99, 'F': 18.47,
                  'Na': 25.69, 'Mg': 34.87, 'Al': 78.23, 'Si': 106.6, 'P': 75.42, 'S': 65.66, 'Cl': 28.59, 'I':24.04}
        # Thermal contribution to enthalpy Hss(298 K) - Hss(0 K) reported by Gaussian thermo whitepaper
        # This will be subtracted from the corresponding value in atomHf to produce an enthalpy used in calculating the enthalpy of formation at 298 K
        atomThermal = {'H': 1.01,
                       'Li': 1.1, 'Be': 0.46, 'B': 0.29, 'C': 0.25, 'N': 1.04, 'O': 1.04, 'F': 1.05,
                       'Na': 1.54, 'Mg': 1.19, 'Al': 1.08, 'Si': 0.76, 'P': 1.28, 'S': 1.05, 'Cl': 1.1, 'I':1.48}
        # Total energy correction used to reach gas-phase reference state
        # Note: Spin orbit coupling no longer included in these energies, since some model chemistries include it automatically
        atomEnthalpyCorrections = {element: atomHf[element] - atomThermal[element] for element in atomHf}
        for symbol, count in atoms.items():
            if symbol in atomEnthalpyCorrections:
                E0 += count * atomEnthalpyCorrections[symbol] * 4184.
            else:
                raise Exception('Element "{}" is not supported.'.format(symbol))

    if applyBondEnergyCorrections:
        # Step 3: Bond energy corrections
        #The order of elements in the bond correction label is important and should follow the order specified below:
        #'C', 'N', 'O', 'S', 'P', and 'H'
        #Use ``-``/``=``/``#`` to denote a single/double/triple bond, respectively.
        # For example, ``'C=N'`` is correct while ``'N=C'`` is incorrect
        bondEnergies = {}
        # 'S-H', 'C-S', 'C=S', 'S-S', 'O-S', 'O=S', 'O=S=O' taken from http://hdl.handle.net/1721.1/98155 (both for
        # 'CCSD(T)-F12/cc-pVDZ-F12' and 'CCSD(T)-F12/cc-pVTZ-F12')
        if modelChemistry == 'ccsd(t)-f12/cc-pvdz-f12':
            bondEnergies = { 'C-H': -0.46, 'C-C': -0.68, 'C=C': -1.90, 'C#C': -3.13,
                'O-H': -0.51, 'C-O': -0.23, 'C=O': -0.69, 'O-O': -0.02, 'C-N': -0.67,
                'C=N': -1.46, 'C#N': -2.79, 'N-O': 0.74, 'N_O': -0.23, 'N=O': -0.51,
                'N-H': -0.69, 'N-N': -0.47, 'N=N': -1.54, 'N#N': -2.05, 'S-H': 0.87,
                'C-S': 0.42, 'C=S': 0.51, 'S-S': 0.86, 'O-S': 0.23, 'O=S': -0.53,
                'O=S=O': 1.95, }
        elif modelChemistry == 'ccsd(t)-f12/cc-pvtz-f12':
            bondEnergies = { 'C-H': -0.09, 'C-C': -0.27, 'C=C': -1.03, 'C#C': -1.79,
                'O-H': -0.06, 'C-O': 0.14, 'C=O': -0.19, 'O-O': 0.16, 'C-N': -0.18,
                'C=N': -0.41, 'C#N': -1.41, 'N-O': 0.87, 'N_O': -0.09, 'N=O': -0.23,
                'N-H': -0.01, 'N-N': -0.21, 'N=N': -0.44, 'N#N': -0.76, 'S-H': 0.52,
                'C-S': 0.13, 'C=S': -0.12, 'S-S': 0.30, 'O-S': 0.15, 'O=S': -2.61,
                'O=S=O': 0.27, }
        elif modelChemistry == 'ccsd(t)-f12/cc-pvqz-f12':
            bondEnergies = { 'C-H': -0.08, 'C-C': -0.26, 'C=C': -1.01, 'C#C': -1.66,
                'O-H':  0.07, 'C-O': 0.25, 'C=O': -0.03, 'O-O': 0.26, 'C-N': -0.20,
                'C=N': -0.30, 'C#N': -1.33, 'N-O': 1.01, 'N_O': -0.03, 'N=O': -0.26,
                'N-H':  0.06, 'N-N': -0.23, 'N=N': -0.37, 'N#N': -0.64,}
        elif modelChemistry == 'cbs-qb3':
            bondEnergies = {
                'C-C': -0.495,'C-H': -0.045,'C=C': -0.825,'C-O': 0.378,'C=O': 0.743,'O-H': -0.423,  #Table2: Paraskevas, PD (2013). Chemistry-A European J., DOI: 10.1002/chem.201301381
                'C#C': -0.64, 'C#N': -0.89, 'C-S': 0.43,  'O=S': -0.78,'S-H': 0.0,  'C-N': -0.13, 'C-Cl': 1.29, 'C-F': 0.55,  # Table IX: Petersson GA (1998) J. of Chemical Physics, DOI: 10.1063/1.477794
                'N-H': -0.42, 'N=O': 1.11,  'N-N': -1.87, 'N=N': -1.58,'N-O': 0.35,  #Table 2: Ashcraft R (2007) J. Phys. Chem. B; DOI: 10.1021/jp073539t
                'N#N': -2.0,  'O=O': -0.2,  'H-H': 1.1,  # Unknown source
                 }
        elif modelChemistry in ['b3lyp/cbsb7', 'b3lyp/6-311g(2d,d,p)', 'b3lyp/6-311+g(3df,2p)', 'b3lyp/6-31g**']:
            bondEnergies = { 'C-H': 0.25, 'C-C': -1.89, 'C=C': -0.40, 'C#C': -1.50,
                'O-H': -1.09, 'C-O': -1.18, 'C=O': -0.01, 'N-H': 1.36, 'C-N': -0.44,
                'C#N': 0.22, 'C-S': -2.35, 'O=S': -5.19, 'S-H': -0.52, }
        else:
            logging.warning('No bond energy correction found for model chemistry: {0}'.format(modelChemistry))

        for symbol, count in bonds.items():
            if symbol in bondEnergies:
                E0 += count * bondEnergies[symbol] * 4184.
            elif symbol[::-1] in bondEnergies:
                E0 += count * bondEnergies[symbol[::-1]] * 4184.
            else:
                logging.warning('Ignored unknown bond type {0!r}.'.format(symbol))

    return E0

class Log(object):
    """
    Represent a general log file.
    The attribute `path` refers to the location on disk of the log file of interest.
    A method is provided to determine whether it is a Gaussian, Molpro, or QChem type.
    """

    def __init__(self, path):
        self.path = path

    def determine_qm_software(self, fullpath):
        """
        Given a path to the log file of a QM software, determine whether it is Gaussian, Molpro, or QChem
        """
        f = open(fullpath, 'r')
        line = f.readline()
        software_log = None
        while line != '':
            if 'gaussian' in line.lower():
                f.close()
                software_log = GaussianLog(fullpath)
                break
            elif 'qchem' in line.lower():
                f.close()
                software_log = QchemLog(fullpath)
                break
            elif 'molpro' in line.lower():
                f.close()
                software_log = MolproLog(fullpath)
                break
            line = f.readline()
        f.close()
        self.software_log = software_log


def projectRotors(conformer, F, rotors, linear, TS):
    """
    For a given `conformer` with associated force constant matrix `F`, lists of
    rotor information `rotors`, `pivots`, and `top1`, and the linearity of the
    molecule `linear`, project out the nonvibrational modes from the force
    constant matrix and use this to determine the vibrational frequencies. The
    list of vibrational frequencies is returned in cm^-1.

    Refer to Gaussian whitepaper (http://gaussian.com/vib/) for procedure to calculate 
    harmonic oscillator vibrational frequencies using the force constant matrix.
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
        with numpy.warnings.catch_warnings():
            numpy.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(numpy.sqrt(eig[i])/(2 * math.pi * constants.c * 100))

    # Now we can start thinking about projecting out the internal rotations
    Dint=numpy.zeros((3*Natoms,Nrotors), numpy.float64)

    counter=0
    for i, rotor in enumerate(rotors):
        if len(rotor) == 5:
            scanLog, pivots, top, symmetry, fit = rotor
        elif len(rotor) == 3:
            pivots, top, symmetry = rotor
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
        with numpy.warnings.catch_warnings():
            numpy.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(numpy.sqrt(eig[i])/(2 * math.pi * constants.c * 100))
        
    return numpy.sqrt(eig[-Nvib:]) / (2 * math.pi * constants.c * 100)

def assign_frequency_scale_factor(model_chemistry):
    """
    Assign the frequency scaling factor according to the model chemistry.
    Refer to https://comp.chem.umn.edu/freqscale/index.html for future updates of these factors
    """
    freq_dict = {'cbs-qb3': 0.99,  # J. Chem. Phys. 1999, 110, 2822–2827
                 # 'g3': ,
                 'm08so/mg3s*': 0.983,  # DOI: 10.1021/ct100326h, taken as 'M08-SO/MG3S'
                 'm06-2x/cc-pvtz': 0.955,  # http://cccbdb.nist.gov/vibscalejust.asp
                 # 'klip_1': ,
                 # 'klip_2': ,
                 # 'klip_3': ,
                 # 'klip_2_cc': ,
                 # 'ccsd(t)-f12/cc-pvdz-f12_h-tz': ,
                 # 'ccsd(t)-f12/cc-pvdz-f12_h-qz': ,
                 'ccsd(t)-f12/cc-pvdz-f12': 0.979,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/cc-pvdz'
                 'ccsd(t)-f12/cc-pvtz-f12': 0.984,  # DOI: 10.1021/ct100326h, taken as 'CCSD(T)-F12a/cc-pVTZ-F12'
                 'ccsd(t)-f12/cc-pvqz-f12': 0.970,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/cc-pvqz'
                 'ccsd(t)-f12/cc-pcvdz-f12': 0.971,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/cc-pcvdz'
                 'ccsd(t)-f12/cc-pcvtz-f12': 0.966,
                 # 'ccsd(t)-f12/cc-pcvqz-f12': ,
                 # 'ccsd(t)-f12/cc-pvtz-f12(-pp)': ,
                 # 'ccsd(t)/aug-cc-pvtz(-pp)': ,
                 'ccsd(t)-f12/aug-cc-pvdz': 0.963,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/aug-cc-pvdz'
                 'ccsd(t)-f12/aug-cc-pvtz': 0.970,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/aug-cc-pvtz'
                 'ccsd(t)-f12/aug-cc-pvqz': 0.975,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as 'ccsd(t)/aug-cc-pvqz'
                 # 'b-ccsd(t)-f12/cc-pvdz-f12': ,
                 # 'b-ccsd(t)-f12/cc-pvtz-f12': ,
                 # 'b-ccsd(t)-f12/cc-pvqz-f12': ,
                 # 'b-ccsd(t)-f12/cc-pcvdz-f12': ,
                 # 'b-ccsd(t)-f12/cc-pcvtz-f12': ,
                 # 'b-ccsd(t)-f12/cc-pcvqz-f12': ,
                 # 'b-ccsd(t)-f12/aug-cc-pvdz': ,
                 # 'b-ccsd(t)-f12/aug-cc-pvtz': ,
                 # 'b-ccsd(t)-f12/aug-cc-pvqz': ,
                 'mp2_rmp2_pvdz': 0.953,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as ',p2/cc-pvdz'
                 'mp2_rmp2_pvtz': 0.950,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as ',p2/cc-pvdz'
                 'mp2_rmp2_pvqz': 0.962,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as ',p2/cc-pvdz'
                 'ccsd-f12/cc-pvdz-f12': 0.947,  # http://cccbdb.nist.gov/vibscalejust.asp, taken as ccsd/cc-pvdz
                 # 'ccsd(t)-f12/cc-pvdz-f12_noscale': ,
                 # 'g03_pbepbe_6-311++g_d_p': ,
                 # 'fci/cc-pvdz': ,
                 # 'fci/cc-pvtz': ,
                 # 'fci/cc-pvqz': ,
                 # 'bmk/cbsb7': ,
                 # 'bmk/6-311g(2d,d,p)': ,
                 'b3lyp/6-31g**': 0.961,  # http://cccbdb.nist.gov/vibscalejust.asp
                 'b3lyp/6-311+g(3df,2p)': 0.967,  # http://cccbdb.nist.gov/vibscalejust.asp
                 }
    scale_factor = freq_dict.get(model_chemistry.lower(), 1)
    if scale_factor == 1:
        logging.warning('No frequency scale factor found for model chemistry {0}; assuming a value of unity.'.format(
            model_chemistry))
    else:
        logging.info('Assigned a frequency scale factor of {0} for model chemistry {1}'.format(
            scale_factor,model_chemistry))
    return scale_factor
