#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

"""
This module provides the :class:`StatMechJob` class, which represents a
statistical mechanics job used to compute and save the statistical mechanics
information for a single species or transition state.
"""

import os.path
import math
import numpy as np
import logging
import subprocess
import os

from rdkit.Chem import GetPeriodicTable

import rmgpy.constants as constants
from rmgpy.species import TransitionState, Species
from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor, FreeRotor
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech import schrodinger
from rmgpy.statmech.mode import Mode
from rmgpy.exceptions import InputError
from rmgpy.quantity import Quantity
from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.element import elementList

from arkane.output import prettify
from arkane.log import Log
from arkane.gaussian import GaussianLog
from arkane.molpro import MolproLog
from arkane.qchem import QChemLog
from arkane.common import symbol_by_number
from arkane.common import ArkaneSpecies

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

        angles = np.array(angles)
        energies = np.array(energies)
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


def hinderedRotor(scanLog, pivots, top, symmetry=None, fit='best'):
    return [scanLog, pivots, top, symmetry, fit]


def freeRotor(pivots,top,symmetry):
    return [pivots,top,symmetry]

def hinderedRotor2D(scandir,pivots1,top1,symmetry1,pivots2,top2,symmetry2,symmetry='none'):
    return [scandir,pivots1,top1,symmetry1,pivots2,top2,symmetry2,symmetry]

class StatMechJob(object):
    """
    A representation of a Arkane statistical mechanics job. This job is used
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
        self.supporting_info = [self.species.label]
        self.bonds = None
        self.arkane_species = ArkaneSpecies(species=species)

    def execute(self, outputFile=None, plot=False, pdep=False):
        """
        Execute the statistical mechanics job, saving the results to the
        given `outputFile` on disk.
        """
        self.load(pdep)
        if outputFile is not None:
            self.save(outputFile)
        logging.debug('Finished statmech job for species {0}.'.format(self.species))
        logging.debug(repr(self.species))

    def load(self, pdep=False):
        """
        Load the statistical mechanics parameters for each conformer from
        the associated files on disk. Creates :class:`Conformer` objects for
        each conformer and appends them to the list of conformers on the
        species object.
        """
        path = self.path
        is_ts = isinstance(self.species, TransitionState)
        _, file_extension = os.path.splitext(path)
        if file_extension in ['.yml', '.yaml']:
            self.arkane_species.load_yaml(path=path, species=self.species, pdep=pdep)
            self.species.conformer = self.arkane_species.conformer
            if is_ts:
                self.species.frequency = self.arkane_species.imaginary_frequency
            else:
                self.species.transportData = self.arkane_species.transport_data
                self.species.energyTransferModel = self.arkane_species.energy_transfer_model
                if self.arkane_species.adjacency_list is not None:
                    self.species.molecule = [Molecule().fromAdjacencyList(adjlist=self.arkane_species.adjacency_list)]
                elif self.arkane_species.inchi is not None:
                    self.species.molecule = [Molecule().fromInChI(inchistr=self.arkane_species.inchi)]
                elif self.arkane_species.smiles is not None:
                    self.species.molecule = [Molecule().fromSMILES(smilesstr=self.arkane_species.smiles)]
            return

        logging.info('Loading statistical mechanics parameters for {0}...'.format(self.species.label))

        global_context = {
            '__builtins__': None,
        }
        local_context = {
            '__builtins__': None,
            'True': True,
            'False': False,
            'HinderedRotor': hinderedRotor,
            'FreeRotor': freeRotor,
            'HinderedRotor2D' : hinderedRotor2D,
            # File formats
            'GaussianLog': GaussianLog,
            'QChemLog': QChemLog,
            'MolproLog': MolproLog,
            'ScanLog': ScanLog,
            'Log': Log
        }

        directory = os.path.abspath(os.path.dirname(path))

        with open(path, 'r') as f:
            try:
                exec f in global_context, local_context
            except (NameError, TypeError, SyntaxError):
                logging.error('The species file {0} was invalid:'.format(path))
                raise

        if self.bonds is None:
            try:
                self.bonds = local_context['bonds']
            except KeyError:
                self.bonds = {}

        try:
            linear = local_context['linear']
        except KeyError:
            logging.error('You did not set whether the molecule is linear with the required `linear` parameter')
            raise

        try:
            externalSymmetry = local_context['externalSymmetry']
        except KeyError:
            externalSymmetry = None

        try:
            spinMultiplicity = local_context['spinMultiplicity']
        except KeyError:
            spinMultiplicity = 0

        try:
            opticalIsomers = local_context['opticalIsomers']
        except KeyError:
            logging.debug('No opticalIsomers provided, estimating them from the quantum file.')
            opticalIsomers = None

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
        E0_withZPE, E0 = None, None
        energyLog = None
        if isinstance(energy, Log) and not isinstance(energy, (GaussianLog,QChemLog,MolproLog)):
            energyLog = determine_qm_software(os.path.join(directory, energy.path))
        elif isinstance(energy, (GaussianLog,QChemLog,MolproLog)):
            energyLog = energy
            energyLog.path = os.path.join(directory, energyLog.path)
        elif isinstance(energy, float):
            E0 = energy
        elif isinstance(energy, tuple) and len(energy) == 2:
            # this is likely meant to be a quantity object with ZPE already accounted for
            energy_temp = Quantity(energy)
            E0_withZPE = energy_temp.value_si # in J/mol
        elif isinstance(energy, tuple) and len(energy) == 3:
            if energy[2] == 'E0':
                energy_temp = Quantity(energy[:2])
                E0 = energy_temp.value_si / constants.E_h / constants.Na# convert J/mol to Hartree
            elif energy[2] == 'E0-ZPE':
                energy_temp = Quantity(energy[:2])
                E0_withZPE = energy_temp.value_si # in J/mol
            else:
                raise InputError('The third argument for E0 energy value should '\
                                 'be E0 (for energy w/o ZPE) or E0-ZPE. Value '\
                                 'entered {0}'.format(energy[2]))
        try:
            geomLog = local_context['geometry']
        except KeyError:
            raise InputError('Required attribute "geometry" not found in species file {0!r}.'.format(path))
        if isinstance(geomLog, Log) and not isinstance(energy, (GaussianLog,QChemLog,MolproLog)):
            geomLog = determine_qm_software(os.path.join(directory, geomLog.path))
        else:
            geomLog.path = os.path.join(directory, geomLog.path)

        try:
            statmechLog = local_context['frequencies']
        except KeyError:
            raise InputError('Required attribute "frequencies" not found in species file {0!r}.'.format(path))
        if isinstance(statmechLog, Log) and not isinstance(energy, (GaussianLog,QChemLog,MolproLog)):
            statmechLog = determine_qm_software(os.path.join(directory, statmechLog.path))
        else:
            statmechLog.path = os.path.join(directory, statmechLog.path)

        if 'frequencyScaleFactor' in local_context:
            logging.warning('Ignoring frequency scale factor in species file {0!r}.'.format(path))

        rotors = []
        if self.includeHinderedRotors:
            try:
                rotors = local_context['rotors']
            except KeyError:
                pass

        # If hindered/free rotors are included in Statmech job, ensure that the same (freq) log file is used for
        # both the species's optimized geometry and Hessian. This approach guarantees that the geometry and Hessian
        # will be defined in the same Cartesian coordinate system ("Input Orientation", as opposed to
        # "Standard Orientation", or something else). Otherwise, if the geometry and Hessian are read from different
        # log files, it is very easy for them to be defined in different coordinate systems, unless the user is very
        # careful. The current implementation only performs this check for Gaussian logs. If QChem logs are used, only
        # a warning is output reminding the user to ensure the geometry and Hessian are defined in consistent
        # coordinates.
        if len(rotors) > 0:
            if isinstance(statmechLog, GaussianLog):
                if statmechLog.path != geomLog.path:
                    raise InputError('For {0!r}, the geometry log, {1!r}, and frequency log, {2!r}, are not the same.'
                                     ' In order to ensure the geometry and Hessian of {0!r} are defined in consistent'
                                     ' coordinate systems for hindered/free rotor projection, either use the frequency'
                                     ' log for both geometry and frequency, or remove rotors.'.format(
                                      self.species.label, geomLog.path, statmechLog.path))
            elif isinstance(statmechLog, QChemLog):
                    logging.warning('QChem log will be used for Hessian of {0!r}. Please verify that the geometry'
                                    ' and Hessian of {0!r} are defined in the same coordinate system'.format(
                                     self.species.label))

        logging.debug('    Reading molecular degrees of freedom...')
        conformer, unscaled_frequencies = statmechLog.loadConformer(symmetry=externalSymmetry,
                                                                    spinMultiplicity=spinMultiplicity,
                                                                    opticalIsomers=opticalIsomers,
                                                                    label=self.species.label)
        for mode in conformer.modes:
            if isinstance(mode, (LinearRotor, NonlinearRotor)):
                self.supporting_info.append(mode)
                break
        if unscaled_frequencies:
            self.supporting_info.append(unscaled_frequencies)

        if conformer.spinMultiplicity == 0:
            raise ValueError("Could not read spin multiplicity from log file {0},\n"
                             "please specify the multiplicity in the input file.".format(self.path))

        logging.debug('    Reading optimized geometry...')
        coordinates, number, mass = geomLog.loadGeometry()

        # Infer atoms from geometry
        atoms = {}
        for atom_num in number:
            try:
                symbol = symbol_by_number[atom_num]
            except KeyError:
                raise Exception('Could not recognize element number {0}.'.format(atom_num))
            atoms[symbol] = atoms.get(symbol, 0) + 1

        # Save atoms for use in writing thermo output
        if isinstance(self.species, Species):
            self.species.props['elementCounts'] = atoms

        conformer.coordinates = (coordinates,"angstroms")
        conformer.number = number
        conformer.mass = (mass,"amu")

        logging.debug('    Reading energy...')
        if E0_withZPE is None:
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
                                        self.bonds,
                                        atomEnergies=self.atomEnergies,
                                        applyAtomEnergyCorrections=self.applyAtomEnergyCorrections,
                                        applyBondEnergyCorrections=self.applyBondEnergyCorrections)
            if len(number) > 1:
                ZPE = statmechLog.loadZeroPointEnergy() * self.frequencyScaleFactor
            else:
                # Monoatomic species don't have frequencies
                ZPE = 0.0
            logging.debug('Corrected minimum energy is {0} J/mol'.format(E0))
            # The E0_withZPE at this stage contains the ZPE
            E0_withZPE = E0 + ZPE

            logging.debug('         Scaling factor used = {0:g}'.format(self.frequencyScaleFactor))
            logging.debug('         ZPE (0 K) = {0:g} kcal/mol'.format(ZPE / 4184.))
            logging.debug('         E0 (0 K) = {0:g} kcal/mol'.format(E0_withZPE / 4184.))

        conformer.E0 = (E0_withZPE*0.001,"kJ/mol")

        # If loading a transition state, also read the imaginary frequency
        if is_ts:
            neg_freq = statmechLog.loadNegativeFrequency()
            self.species.frequency = (neg_freq * self.frequencyScaleFactor, "cm^-1")
            self.supporting_info.append(neg_freq)

        # Read and fit the 1D hindered rotors if applicable
        # If rotors are found, the vibrational frequencies are also
        # recomputed with the torsional modes removed

        F = statmechLog.loadForceConstantMatrix()

        if F is not None and len(mass) > 1 and len(rotors) > 0:

            logging.debug('    Fitting {0} hindered rotors...'.format(len(rotors)))
            rotorCount = 0
            for j,q in enumerate(rotors):
                symmetry = None
                if len(q) == 3:
                    # No potential scan is given, this is a free rotor
                    pivots, top, symmetry = q
                    inertia = conformer.getInternalReducedMomentOfInertia(pivots, top) * constants.Na * 1e23
                    rotor = FreeRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                    conformer.modes.append(rotor)
                    rotorCount += 1
                elif len(q) == 8:
                    scandir,pivots1,top1,symmetry1,pivots2,top2,symmetry2,symmetry = q
                    logging.info("Calculating energy levels for 2D-HR, may take a while...")
                    rotor = HinderedRotor2D(name='r'+str(j),torsigma1=symmetry1,torsigma2=symmetry2,symmetry=symmetry,
                                            calcPath=os.path.join(directory,scandir),pivots1=pivots1,pivots2=pivots2,top1=top1,top2=top2)
                    rotor.run()
                    conformer.modes.append(rotor)
                    rotorCount += 2
                elif len(q) in [4, 5]:
                    # This is a hindered rotor
                    if len(q) == 5:
                        scanLog, pivots, top, symmetry, fit = q
                    elif len(q) == 4:
                        # the symmetry number will be derived from the scan
                        scanLog, pivots, top, fit = q
                    # Load the hindered rotor scan energies
                    if isinstance(scanLog, Log) and not isinstance(energy, (GaussianLog,QChemLog,MolproLog)):
                        scanLog = determine_qm_software(os.path.join(directory, scanLog.path))
                    if isinstance(scanLog, GaussianLog):
                        scanLog.path = os.path.join(directory, scanLog.path)
                        v_list, angle = scanLog.loadScanEnergies()
                        scanLogOutput = ScanLog(os.path.join(directory, '{0}_rotor_{1}.txt'.format(
                            self.species.label, rotorCount+1)))
                        scanLogOutput.save(angle, v_list)
                    elif isinstance(scanLog, QChemLog):
                        scanLog.path = os.path.join(directory, scanLog.path)
                        v_list, angle = scanLog.loadScanEnergies()
                        scanLogOutput = ScanLog(os.path.join(directory, '{0}_rotor_{1}.txt'.format(
                            self.species.label, rotorCount+1)))
                        scanLogOutput.save(angle, v_list)
                    elif isinstance(scanLog, ScanLog):
                        scanLog.path = os.path.join(directory, scanLog.path)
                        angle, v_list = scanLog.load()
                    else:
                        raise InputError('Invalid log file type {0} for scan log.'.format(scanLog.__class__))

                    if symmetry is None:
                        symmetry = determine_rotor_symmetry(v_list, self.species.label, pivots)
                    inertia = conformer.getInternalReducedMomentOfInertia(pivots, top) * constants.Na * 1e23

                    cosineRotor = HinderedRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                    cosineRotor.fitCosinePotentialToData(angle, v_list)
                    fourierRotor = HinderedRotor(inertia=(inertia,"amu*angstrom^2"), symmetry=symmetry)
                    fourierRotor.fitFourierPotentialToData(angle, v_list)

                    Vlist_cosine = np.zeros_like(angle)
                    Vlist_fourier = np.zeros_like(angle)
                    for i in range(angle.shape[0]):
                        Vlist_cosine[i] = cosineRotor.getPotential(angle[i])
                        Vlist_fourier[i] = fourierRotor.getPotential(angle[i])

                    if fit == 'cosine':
                        rotor = cosineRotor
                        rotorCount += 1
                        conformer.modes.append(rotor)
                    elif fit == 'fourier':
                        rotor = fourierRotor
                        rotorCount += 1
                        conformer.modes.append(rotor)
                    elif fit == 'best':
                        rms_cosine = np.sqrt(np.sum((Vlist_cosine - v_list) * (Vlist_cosine - v_list)) /
                                                (len(v_list) - 1)) / 4184.
                        rms_fourier = np.sqrt(np.sum((Vlist_fourier - v_list) * (Vlist_fourier - v_list))/
                                                 (len(v_list) - 1)) / 4184.

                        # Keep the rotor with the most accurate potential
                        rotor = cosineRotor if rms_cosine < rms_fourier else fourierRotor
                        # However, keep the cosine rotor if it is accurate enough, the
                        # fourier rotor is not significantly more accurate, and the cosine
                        # rotor has the correct symmetry
                        if rms_cosine < 0.05 and rms_cosine / rms_fourier < 2.0 and rms_cosine / rms_fourier < 4.0\
                                and symmetry == cosineRotor.symmetry:
                            rotor = cosineRotor

                        conformer.modes.append(rotor)

                        self.plotHinderedRotor(angle, v_list, cosineRotor, fourierRotor, rotor, rotorCount, directory)

                        rotorCount += 1

            logging.debug('    Determining frequencies from reduced force constant matrix...')
            frequencies = np.array(projectRotors(conformer, F, rotors, linear, is_ts))

        elif len(conformer.modes) > 2:
            if len(rotors) > 0:
                logging.warn('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already'
                             ' present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for'
                             ' Gaussian to generate the force constant matrix, if running Molpro include keyword print,'
                             ' hessian')
            frequencies = conformer.modes[2].frequencies.value_si
            rotors = np.array([])
        else:
            if len(rotors) > 0:
                logging.warn('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already'
                             ' present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for'
                             ' Gaussian to generate the force constant matrix, if running Molpro include keyword print,'
                             ' hessian')
            frequencies = np.array([])
            rotors = np.array([])

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
            f.write('#   {0} {1:9.4f} {2:9.4f} {3:9.4f}\n'.format(symbol_by_number[number[i]], x, y, z))

        result = 'conformer(label={0!r}, E0={1!r}, modes={2!r}, spinMultiplicity={3:d}, opticalIsomers={4:d}'.format(
            self.species.label,
            conformer.E0,
            conformer.modes,
            conformer.spinMultiplicity,
            conformer.opticalIsomers,
        )
        try:
            result += ', frequency={0!r}'.format(self.species.frequency)
        except AttributeError: pass
        result += ')'
        f.write('{0}\n\n'.format(prettify(result)))
        f.close()

    def plotHinderedRotor(self, angle, v_list, cosineRotor, fourierRotor, rotor, rotorIndex, directory):
        """
        Plot the potential for the rotor, along with its cosine and Fourier
        series potential fits. The plot is saved to a set of files of the form
        ``hindered_rotor_1.pdf``.
        """
        try:
            import pylab
        except ImportError:
            return

        phi = np.arange(0, 6.3, 0.02, np.float64)
        Vlist_cosine = np.zeros_like(phi)
        Vlist_fourier = np.zeros_like(phi)
        for i in range(phi.shape[0]):
            Vlist_cosine[i] = cosineRotor.getPotential(phi[i])
            Vlist_fourier[i] = fourierRotor.getPotential(phi[i])

        fig = pylab.figure(figsize=(6,5))
        pylab.plot(angle, v_list / 4184., 'ok')
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
    H doublet, C triplet, O triplet, N quartet, S triplet, P quartet, I doublet.

    `bonds` is a dictionary associating bond types with the number
    of that bond in the molecule.
    """

    if applyAtomEnergyCorrections:
        # Spin orbit correction (SOC) in Hartrees
        # Values taken from ref 22 of http://dx.doi.org/10.1063/1.477794 and converted to hartrees
        # Values in millihartree are also available (with fewer significant figures) from table VII of http://dx.doi.org/10.1063/1.473182
        # Iodine SOC calculated as a weighted average of the electronic spin splittings of the lowest energy state. The splittings are
        # obtained from Huber, K.P.; Herzberg, G., Molecular Spectra and Molecular Structure. IV. Constants of Diatomic Molecules, Van Nostrand Reinhold Co., 1979
        SOC = {'H': 0.0, 'N': 0.0, 'O': -0.000355, 'C': -0.000135, 'S': -0.000893, 'P': 0.0, 'I':-0.011547226,}

        # Step 1: Reference all energies to a model chemistry-independent basis
        # by subtracting out that model chemistry's atomic energies
        # All model chemistries here should be lower-case because the user input is changed to lower-case
        if atomEnergies is None:
            # Note: If your model chemistry does not include spin orbit coupling, you should add the corrections to the energies here
            if modelChemistry.startswith('cbs-qb3'):  # only check start of string to allow different bond corrections (see below)
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

            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvdz':  # note that all atom corrections but S are fitted, the correction for S is calculated
                atomEnergies = {'H':-0.499459066131 + SOC['H'], 'N':-54.524279516472 + SOC['N'], 'O':-74.992097308083+ SOC['O'], 'C':-37.786694171716+ SOC['C'], 'S':-397.648733842400 + SOC['S']}
            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvtz':
                atomEnergies = {'H':-0.499844820798 + SOC['H'], 'N':-54.527419359906 + SOC['N'], 'O':-75.000001429806 + SOC['O'], 'C':-37.788504810868 + SOC['C'], 'S':-397.666903000231 + SOC['S']}
            elif modelChemistry == 'ccsd(t)-f12/aug-cc-pvqz':
                atomEnergies = {'H':-0.499949526073 + SOC['H'], 'N':-54.529569719016 + SOC['N'], 'O':-75.004026586610+ SOC['O'], 'C':-37.789387892348+ SOC['C'], 'S':-397.671214204994 + SOC['S']}


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
            elif modelChemistry == 'wb97x-d/aug-cc-pvtz':
                atomEnergies = {'H':-0.502803+ SOC['H'], 'N':-54.585652+ SOC['N'], 'O':-75.068286+ SOC['O'], 'C':-37.842014+ SOC['C']}

            elif modelChemistry == 'MRCI+Davidson/aug-cc-pV(T+d)Z':  # Calculated atomic energies (unfitted)
                atomEnergies = {'H':-0.49982118 + SOC['H'], 'C':-37.78321274 + SOC['C'], 'N':-54.51729444 + SOC['N'], 'O':-74.97847534 + SOC['O'], 'S':-397.6571654 + SOC['S']}

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
        # Experimental enthalpy of formation at 0 K, 1 bar for gas phase
        # See Gaussian thermo whitepaper at http://www.gaussian.com/g_whitepap/thermo.htm)
        # Note: These values are relatively old and some improvement may be possible by using newer values
        # (particularly for carbon).
        # However, care should be taken to ensure that they are compatible with the BAC values (if BACs are used)
        # He, Ne, K, Ca, Ti, Cu, Zn, Ge, Br, Kr, Rb, Ag, Cd, Sn, I, Xe, Cs, Hg, and Pb are taken from CODATA
        # Codata: Cox, J. D., Wagman, D. D., and Medvedev, V. A., CODATA Key Values for Thermodynamics, Hemisphere
        # Publishing Corp., New York, 1989. (http://www.science.uwaterloo.ca/~cchieh/cact/tools/thermodata.html)
        atom_hf = {'H': 51.63, 'He': -1.481,
                  'Li': 37.69, 'Be': 76.48, 'B': 136.2, 'C': 169.98, 'N': 112.53, 'O': 58.99, 'F': 18.47, 'Ne': -1.481,
                  'Na': 25.69, 'Mg': 34.87, 'Al': 78.23, 'Si': 106.6, 'P': 75.42, 'S': 65.66, 'Cl': 28.59,
                  'K': 36.841, 'Ca': 41.014, 'Ti': 111.2, 'Cu': 79.16, 'Zn': 29.685, 'Ge': 87.1, 'Br': 25.26, 'Kr': -1.481,
                  'Rb': 17.86, 'Ag': 66.61, 'Cd': 25.240, 'Sn': 70.50, 'I': 24.04, 'Xe': -1.481,
                  'Cs': 16.80, 'Hg': 13.19, 'Pb': 15.17}
        # Thermal contribution to enthalpy Hss(298 K) - Hss(0 K) reported by Gaussian thermo whitepaper
        # This will be subtracted from the corresponding value in atom_hf to produce an enthalpy used in calculating
        # the enthalpy of formation at 298 K
        atom_thermal = {'H': 1.01, 'He': 1.481,
                       'Li': 1.1, 'Be': 0.46, 'B': 0.29, 'C': 0.25, 'N': 1.04, 'O': 1.04, 'F': 1.05, 'Ne': 1.481,
                       'Na': 1.54, 'Mg': 1.19, 'Al': 1.08, 'Si': 0.76, 'P': 1.28, 'S': 1.05, 'Cl': 1.1,
                       'K': 1.481, 'Ca': 1.481, 'Ti': 1.802, 'Cu': 1.481, 'Zn': 1.481, 'Ge': 1.768, 'Br': 1.481, 'Kr': 1.481,
                       'Rb': 1.481, 'Ag': 1.481, 'Cd': 1.481, 'Sn': 1.485, 'I': 1.481, 'Xe': 1.481,
                       'Cs': 1.481, 'Hg': 1.481, 'Pb': 1.481}
        # Total energy correction used to reach gas-phase reference state
        # Note: Spin orbit coupling is no longer included in these energies, since some model chemistries include it
        # automatically
        atom_enthalpy_corrections = {element: atom_hf[element] - atom_thermal[element] for element in atom_hf}
        for symbol, count in atoms.items():
            if symbol in atom_enthalpy_corrections:
                E0 += count * atom_enthalpy_corrections[symbol] * 4184.
            else:
                raise Exception('Element "{0}" is not yet supported in Arkane.'
                                ' To include it, add its experimental heat of formation'.format(symbol))

    if applyBondEnergyCorrections:
        # Step 3: Bond energy corrections
        # The order of elements in the bond correction label is important and should follow the order specified below:
        # 'C', 'N', 'O', 'S', 'P', and 'H'
        # Use ``-``/``=``/``#`` to denote a single/double/triple bond, respectively.
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
                'C-H': -0.11, 'C-C': -0.30, 'C=C': -0.08, 'C#C': -0.64, 'O-H' : 0.02, 'C-O': 0.33, 'C=O': 0.55,  # Table IX: Petersson GA (1998) J. of Chemical Physics, DOI: 10.1063/1.477794
                'N-H': -0.42, 'C-N': -0.13, 'C#N': -0.89, 'C-F':  0.55, 'C-Cl': 1.29, 'S-H': 0.0,  'C-S': 0.43, 'O=S': -0.78,
                'N=O':  1.11, 'N-N': -1.87, 'N=N': -1.58, 'N-O':  0.35,  #Table 2: Ashcraft R (2007) J. Phys. Chem. B; DOI: 10.1021/jp073539t
                'N#N':  -2.0, 'O=O': -0.2,  'H-H': 1.1,  # Unknown source
            }
        elif modelChemistry == 'cbs-qb3-paraskevas':  # NOTE: The Paraskevas corrections are inaccurate for non-oxygenated hydrocarbons, and may do poorly in combination with the Petersson corrections
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

def determine_qm_software(fullpath):
    """
    Given a path to the log file of a QM software, determine whether it is Gaussian, Molpro, or QChem
    """
    with open(fullpath, 'r') as f:
        line = f.readline()
        software_log = None
        while line != '':
            if 'gaussian' in line.lower():
                f.close()
                software_log = GaussianLog(fullpath)
                break
            elif 'qchem' in line.lower():
                f.close()
                software_log = QChemLog(fullpath)
                break
            elif 'molpro' in line.lower():
                f.close()
                software_log = MolproLog(fullpath)
                break
            line = f.readline()
        else:
            raise InputError("File at {0} could not be identified as a Gaussian, QChem or Molpro log file.".format(fullpath))
    return software_log


def projectRotors(conformer, F, rotors, linear, is_ts):
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
    Nvib = 3 * Natoms - (5 if linear else 6) - Nrotors - (1 if (is_ts) else 0)
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
    amass=np.sqrt(mass/constants.amu)

    # Rotation matrix
    I=conformer.getMomentOfInertiaTensor()
    PMoI, Ixyz = np.linalg.eigh(I)

    external=6
    if linear:
        external=5

    D = np.zeros((Natoms*3,external), np.float64)

    P = np.zeros((Natoms,3), np.float64)

    # Transform the coordinates to the principal axes
    P = np.dot(coordinates,Ixyz)

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

    I = np.identity(Natoms*3, np.float64)

    P = np.zeros((Natoms*3,3*Natoms+external), np.float64)

    P[:,0:external] = D[:,0:external]
    P[:,external:external+3*Natoms] = I[:,0:3*Natoms]

    for i in range(3*Natoms+external):
        norm=0.0
        for j in range(3*Natoms):
            norm+=P[j,i]*P[j,i]
        for j in range(3*Natoms):
            if (norm>1E-15):
                P[j,i]/=np.sqrt(norm)
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
    T = np.zeros((Natoms*3,3*Natoms-external), np.float64)

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

    Fint = np.dot(T.T, np.dot(Fm,T))

    # Get eigenvalues of internal force constant matrix, V = 3N-6 * 3N-6
    eig, V = np.linalg.eigh(Fint)

    logging.debug('Frequencies from internal Hessian')
    for i in range(3*Natoms-external):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i])/(2 * math.pi * constants.c * 100))

    # Now we can start thinking about projecting out the internal rotations
    Dint=np.zeros((3*Natoms,Nrotors), np.float64)

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
    Vmw = np.dot(T,V)
    eigM = np.zeros((3*Natoms-external,3*Natoms-external), np.float64)

    for i in range(3*Natoms-external):
        eigM[i,i]=eig[i]

    Fm=np.dot(Vmw,np.dot(eigM,Vmw.T))

    # Internal rotations are not normal modes => project them on the normal modes and orthogonalize
    # Dintproj =  (3N-6) x (3N) x (3N) x (Nrotors)
    Dintproj=np.dot(Vmw.T,Dint)

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
            Dint[j,i]/=np.sqrt(norm)
        for j in range(i+1,Nrotors):
            proj=0.0
            for k in range (3*Natoms):
                proj+=Dint[k,i]*Dint[k,j]
            for k in range(3*Natoms):
                Dint[k,j]-=proj*Dint[k,i]

    Dintproj=np.dot(Vmw.T,Dint)
    Proj = np.dot(Dint, Dint.T)
    I = np.identity(Natoms*3, np.float64)
    Proj = I - Proj
    Fm=np.dot(Proj, np.dot(Fm,Proj))
    # Get eigenvalues of mass-weighted force constant matrix
    eig, V = np.linalg.eigh(Fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation

    logging.debug('Frequencies from projected Hessian')
    for i in range(3*Natoms):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i])/(2 * math.pi * constants.c * 100))

    return np.sqrt(eig[-Nvib:]) / (2 * math.pi * constants.c * 100)


def assign_frequency_scale_factor(model_chemistry):
    """
    Assign the frequency scaling factor according to the model chemistry.
    Refer to https://comp.chem.umn.edu/freqscale/index.html for future updates of these factors
    """
    freq_dict = {'cbs-qb3': 0.99,  # J. Chem. Phys. 1999, 110, 2822–2827
                 'cbs-qb3-paraskevas': 0.99,
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
                 'ccsd(t)-f12/cc-pvtz-f12': 0.984,  # Taken from https://comp.chem.umn.edu/freqscale/version3b2.htm as CCSD(T)-F12a/cc-pVTZ-F12
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
                 'wb97x-d/aug-cc-pvtz': 0.974, # Taken from https://comp.chem.umn.edu/freqscale/version3b2.htm as ωB97X-D/maug-cc-pVTZ
                 }
    scale_factor = freq_dict.get(model_chemistry.lower(), 1)
    if scale_factor == 1:
        logging.warning('No frequency scale factor found for model chemistry {0}; assuming a value of unity.'.format(
            model_chemistry))
    else:
        logging.info('Assigned a frequency scale factor of {0} for model chemistry {1}'.format(
            scale_factor,model_chemistry))
    return scale_factor


def determine_rotor_symmetry(energies, label, pivots):
    """
    Determine the rotor symmetry number from the potential scan given in :list:`energies` in J/mol units
    Assumes the list represents a 360 degree scan
    str:`label` is the species name, used for logging and error messages
    list:`pivots` are the rotor's pivots, used for logging and error messages
    The *worst* resolution for each peak and valley is determined.
    The first criterion for a symmetric rotor is that the highest peak and the lowest peak must be within the
    worst peak resolution (and the same is checked for valleys).
    A second criterion for a symmetric rotor is that the highest and lowest peaks must be within 10% of
    the highest peak value. This is only applied if the highest peak is above 2 kJ/mol.
    """
    symmetry = None
    min_e = min(energies)
    max_e = max(energies)
    if max_e > 2000:
        tol = 0.10 * max_e  # tolerance for the second criterion
    else:
        tol = max_e
    peaks, valleys = list(), [energies[0]]  # the peaks and valleys of the scan
    worst_peak_resolution, worst_valley_resolution = 0, max(energies[1] - energies[0], energies[-2] - energies[-1])
    for i, e in enumerate(energies):
        # identify peaks and valleys, and determine worst resolutions in the scan
        if i != 0 and i != len(energies) - 1:
            last_point = energies[i - 1]
            next_point = energies[i + 1]
            # this is an intermediate point in the scan
            if e > last_point and e > next_point:
                # this is a local peak
                if any([diff > worst_peak_resolution for diff in [e - last_point, e - next_point]]):
                    worst_peak_resolution = max(e - last_point, e - next_point)
                peaks.append(e)
            elif e < last_point and e < next_point:
                # this is a local valley
                if any([diff > worst_valley_resolution for diff in [energies[i - 1] - e, next_point - e]]):
                    worst_valley_resolution = max(last_point - e, next_point - e)
                valleys.append(e)
    # The number of peaks and valley must always be the same (what goes up must come down), if it isn't then there's
    # something seriously wrong with the scan
    if len(peaks) != len(valleys):
        raise InputError('Rotor of species {0} between pivots {1} does not have the same number'
                         ' of peaks and valleys.'.format(label, pivots))
    min_peak = min(peaks)
    max_peak = max(peaks)
    min_valley = min(valleys)
    max_valley = max(valleys)
    # Criterion 1: worst resolution
    if max_peak - min_peak > worst_peak_resolution:
        # The rotor cannot be symmetric
        symmetry = 1
        reason = 'worst peak resolution criterion'
    elif max_valley - min_valley > worst_valley_resolution:
        # The rotor cannot be symmetric
        symmetry = 1
        reason = 'worst valley resolution criterion'
    # Criterion 2: 10% * max_peak
    elif max_peak - min_peak > tol:
        # The rotor cannot be symmetric
        symmetry = 1
        reason = '10% of the maximum peak criterion'
    else:
        # We declare this rotor as symmetric and the symmetry number in the number of peaks (and valleys)
        symmetry = len(peaks)
        reason = 'number of peaks and valleys, all within the determined resolution criteria'
    if symmetry not in [1, 2, 3]:
        logging.warn('Determined symmetry number {0} for rotor of species {1} between pivots {2};'
                     ' you should make sure this makes sense'.format(symmetry, label, pivots))
    else:
        logging.info('Determined a symmetry number of {0} for rotor of species {1} between pivots {2}'
                     ' based on the {3}.'.format(symmetry, label, pivots, reason))
    return symmetry

class HinderedRotor2D(Mode):
   """
   A statistical mechanical model of a 2D-dimensional hindered rotor.
   Computation of eigenvalues and fourier fitting outsourced to
   Q2DTor software

   The attributes are:

   ======================== ===================================================
   Attribute                Description
   ======================== ===================================================
   `name`                   The Q2DTor name of the rotor
   `torsion1`               The 1-indexed atom indices of the atoms involved in the first rotor
   `torsion2`               The 1-indexed atom indices of the atoms involved in the second rotor
   `torsigma1`              The symmetry number of the first rotor
   `torsigma2`              The symmetry number of the second rotor
   `calcPath`               The directory containing all of the rotor point calculations formated:  name_angle1_angle2
   `symmetry`               Q2DTor symmetry identifier ('none','a','b','c','ab','bc','ac','abc')
   `evals`                  Array of energy levels for the 2D-HR
   `energy`                 Function mapping quantum number to energy level
   ======================== ===================================================
   """

   def __init__(self,name,torsigma1,torsigma2,calcPath,symmetry='none',pivots1=None,pivots2=None,top1=None,top2=None):
       Mode.__init__(self, True)
       self.dof = 2

       self.name = name
       self.calcPath = calcPath

       self.torsion1 = None
       self.torsion2 = None
       self.torsigma1 = torsigma1
       self.torsigma2 = torsigma2
       self.symmetry = symmetry

       self.pivots1 = pivots1
       self.pivots2 = pivots2
       self.top1 = top1
       self.top2 = top2

       self.xyzs = []
       self.phi1s = []
       self.phi2s = []
       self.Es = []
       self.atnums = []

       #prepare a directory to run q2dtor in
       self.q2dtor_path = os.path.join(os.path.split(calcPath)[0],name)
       try:
           os.mkdir(self.q2dtor_path)
           os.mkdir(os.path.join(self.q2dtor_path,'IOfiles'))
       except OSError:
           pass

   def getTorsions(self):
       """
       determine torsions, not entirely necessary for 2D-NS, but important for E2DT
       """
       if not self.torsion1:
           self.readGjf() #check it there is a gaussian format file

       Natoms = len(self.xyzs[0]) #define a feasible torsion from pivots and tops
       aset = set(list(xrange(1,Natoms+1)))
       if not self.torsion1 and self.pivots1 and self.top1:
           if self.pivots1[0] in self.top1:
               self.torsion1 = [list(aset-set(self.top1))[0],self.pivots1[0],self.pivots1[1],self.top1[0]]
           else:
               self.torsion1 = [list(aset-set(self.top1))[0],self.pivots1[1],self.pivots1[0],self.top1[0]]

       if not self.torsion2 and self.pivots2 and self.top2:
           if self.pivots2[0] in self.top2:
               self.torsion2 = [list(aset-set(self.top2))[0],self.pivots2[0],self.pivots2[1],self.top2[0]]
           else:
               self.torsion2 = [list(aset-set(self.top2))[0],self.pivots2[1],self.pivots2[0],self.top2[0]]

   def readScan(self):
       """
       Read quantum optimization job files at self.calcPath to determine
       vectors of angles (self.phi1s, self.phi2s), xyz coordinates (self.xyzs)
       energies (self.Es) and atom numbers (self.atnums) for each point
       """
       phi1s = []
       phi2s = []
       xyzs = []
       Es = []
       atnums = []
       for f in os.listdir(self.calcPath):
           if len(f.split('_')) != 4:
               continue
           s,name,phi1,phi2 = f.split('_')
           phi2,idf = '.'.join(phi2.split('.')[:-1]),phi2.split('.')[-1]
           if idf != 'out':
               continue
           phi1s.append(float(phi1))
           phi2s.append(float(phi2.split(".")[0]))

           fpath = os.path.join(self.calcPath,f)
           lg = Log(fpath)
           lg.determine_qm_software(fpath)

           Es.append(lg.software_log.loadEnergy())
           xyz,atnum,_ = lg.software_log.loadGeometry()
           xyzs.append(xyz)
           atnums.append(atnum)

       self.xyzs = xyzs
       self.phi1s = phi1s
       self.phi2s = phi2s
       self.Es = Es
       self.atnums = atnums

   def readGjf(self):
       """
       read gaussian input file to determine torsions, charge and multiplicity
       unnecessary for 2D-NS
       """
       self.torsion1 = None
       self.torsion2 = None
       self.charge = None
       self.multiplicity = None

       for f in os.listdir(self.calcPath):
           if len(f.split('_')) != 4:
               continue
           s,name,phi1,phi2 = f.split('_')
           phi2,idf = phi2.split('.')
           if idf == 'gjf' and float(phi1) == 0.0 and float(phi2) == 0.0:
               fop = open(os.path.join(self.calcPath,f),'rU')
               lines = fop.readlines()
               for i,line in enumerate(lines):
                   splitLine = line.split()
                   if len(splitLine) < 2:
                       continue
                   elif splitLine[-1] == 'F' and len(splitLine) == 5:
                       if not self.torsion1:
                           self.torsion1 = [int(splitLine[i]) for i in xrange(4)]
                       elif not self.torsion2:
                           self.torsion2 = [int(splitLine[i]) for i in xrange(4)]
                   elif not self.charge and len(splitLine) == 2 and len(lines[i+1].split()) == 4 and len(lines[i+1].split()[0]) <= 2:
                       self.charge = int(splitLine[0])
                       self.multiplicity = int(splitLine[1])



   def writeXYZ(self):
       """
       write an .xyz file for Q2DTor
       done based on the angle coordinates (0.0,0.0)
       """
       atdict = {x.number:x.symbol for x in elementList}
       for i in xrange(len(self.phi1s)):
           if self.phi1s[i] == 0.0 and self.phi2s[i] == 0.0:
               f = open(os.path.join(self.q2dtor_path,self.name+".xyz"),'w')
               f.write(str(len(self.atnums[i]))+'\n')
               f.write("reference geometry for {0}\n".format(self.name))
               element_names = [atdict[k] for k in self.atnums[i]]
               for j,ename in enumerate(element_names):
                   f.write('{ename}    {x}   {y}   {z}\n'.format(ename=ename,
                           x=self.xyzs[i][j][0],y=self.xyzs[i][j][1],z=self.xyzs[i][j][2]))
               f.close()
               break

   def writePes(self):
       """
       write a .pes file for Q2DTor based on the
       read in scans
       """
       atdict = {x.number:x.symbol for x in elementList}
       assert len(self.Es) > 0
       f = open(os.path.join(self.q2dtor_path,'IOfiles',self.name+".pes"),'w')
       for i in xrange(len(self.phi1s)):
           f.write(str(len(self.atnums[i]))+'\n')
           f.write("Geometry   {E}   {phi1}   {phi2}  {name}_{phi1}_{phi2}  YES\n".format(E=self.Es[i]/(constants.E_h*constants.Na),
                   phi1=self.phi1s[i],phi2=self.phi2s[i],name=self.name))
           element_names = [atdict[k] for k in self.atnums[i]]
           for j,ename in enumerate(element_names):
               f.write('{ename}    {x}   {y}   {z}\n'.format(ename=ename,
                           x=self.xyzs[i][j][0],y=self.xyzs[i][j][1],z=self.xyzs[i][j][2]))
       f.close()

   def writeInp(self):
       """
       Write an input file for Q2DTor based on object
       information
       """
       inp = """#----------------------------------#
# Torsional information            #
#----------------------------------#
start_torsions
   torsion1         {tor1}         # atoms involved in torsion 1
   torsion2         {tor2}         # atoms involved in torsion 2
   tsigma1          {torsigma1}               # torsional symmetry number of hindered rotor 1
   tsigma2          {torsigma2}               # torsional symmetry number of hindered rotor 2
end_torsions
#----------------------------------#
# Calculations                     #
#----------------------------------#
start_calcs
   level            hf 3-21g        # the calculation level
   charge           0               # charge
   multiplicity     0               # spin multiplicity, shouldn't matter in context of Arkane
end_calcs
#----------------------------------#
# Torsional PES                    #
#----------------------------------#
start_pes
   t1step           10.0            # step in phi1 for scan calculation [degrees]
   t2step           10.0            # step in phi2 for scan calculation [degrees]
   symmetry         {symmetry}               # Symmetry condition for PES: [a,b,c,ab,ac,bc,abc] or none
end_pes
#----------------------------------#
# Fitting details                  #
#----------------------------------#
start_fourier
   weight           0.9             #
   ignore           0.0             # Set to zero coefficients with smaller absolute value (in cm^-1)
   # Fourier Terms (Even)           #
   cos1             1-9             # i values in cos(i*Phi_1)
   cos2             1-9             # j values in cos(j*Phi_2)
   cos1cos2         1-7 , 1-7       # i,j values in cos(i*Phi_1) * cos(j*Phi_2)
   sin1sin2         1-7 , 1-7       # i,j values in sin(i*Phi_1) * sin(j*Phi_2)
   # Fourier Terms (Odd)            #
   sin1             1-9            # i values in sin(i*Phi_1)
   sin2             1-9            # j values in sin(j*Phi_2)
   cos1sin2         1-7 , 1-7            # i,j values in cos(i*Phi_1) * sin(j*Phi_2)
   sin1cos2         1-7 , 1-7            # i,j values in sin(i*Phi_1) * cos(j*Phi_2)
end_fourier
#----------------------------------#
# Search and Opt stationary points #
#----------------------------------#
start_statpoint
   tolerance        1.0             # step (in degrees) to explore torsional PES when looking for CPs
   freqscal         1.000           # scaling factor for frequencies
end_statpoint
#----------------------------------#
# 2D-NS Hamiltonian                #
#----------------------------------#
start_tor2dns
   dijvar           yes             # yes (dij not constant) or no (dij constant)
   kmax             100             # check 2013-JChemPhys_138_134112, eq (14)
   maxeigen         1e4             # threshold for H eigenvalues (in cm^-1)
end_tor2dns
#----------------------------------#
# Partition functions              #
#----------------------------------#
start_rovibpf
   interpolation    fourier         # fourier or spline order (1,3,5)
   integrationstep  1.0             # integration dphi
end_rovibpf
#----------------------------------#
# Working temperatures             #
#----------------------------------#
start_temperatures                 #
   100.0   150.0   200.0          #
   250.0   300.0   400.0          #
   500.0   700.0  1000.0          #
   1500.0  2000.0  2500.0         #
end_temperatures                   #
#----------------------------------#""".format(tor1="-".join([str(x) for x in self.torsion1]),
                  tor2="-".join([str(x) for x in self.torsion2]),
           torsigma1=self.torsigma1,torsigma2=self.torsigma2,
           charge=self.charge,multiplicity=self.multiplicity,
           symmetry=self.symmetry)
       f = open(os.path.join(self.q2dtor_path,self.name+".inp"),'w')
       f.write(inp)
       f.close()

   def getIcsFile(self):
       """
       use Q2DTor to generate a .ics file
       """
       out = subprocess.check_call(['python2',os.environ['Q2DTor'],self.name,'--init'],
                       cwd=self.q2dtor_path)

   def fitFourier(self):
       """
       use Q2DTor to fit fourier coefficients
       to the potential
       """
       out = subprocess.check_call(['python2',os.environ['Q2DTor'],self.name,'--fourier'],
                       cwd=self.q2dtor_path)

   def getSplistfile(self):
       """
       use Q2DTor to generate a .splist file
       """
       out = subprocess.check_call(['python2',os.environ['Q2DTor'],self.name,'--findsp'],
                       cwd=self.q2dtor_path)

   def getEigvals(self):
       """
       use Q2DTor to determine the QM energy levels for the 2D-NS
       rotors
       writes a .evals file and reads it to fill self.evals and self.energy
       """
       out = subprocess.check_call(['python2',os.environ['Q2DTor'],self.name,'--tor2dns'],
                       cwd=self.q2dtor_path)
       self.readEigvals()

   def readEigvals(self):
       """
       reads an available .evals file to get the QM energy levels
       for the 2D-NS rotors
       """
       f = open(os.path.join(self.q2dtor_path,'IOfiles',self.name+'.evals'),'rU')
       out = f.readlines()
       evals = [float(x.split()[1]) for x in out[2:]] #cm^-1
       self.evals = np.array(evals)*10**2*constants.c*constants.h*constants.Na #J/mol
       self.energy = lambda x: self.evals[x]

   def getPartitionFunction(self,T):
       return schrodinger.getPartitionFunction(T,self.energy,nmax=len(self.evals))

   def getHeatCapacity(self,T):
       return schrodinger.getHeatCapacity(T,self.energy,nmax=len(self.evals))

   def getEnthalpy(self,T):
       return schrodinger.getEnthalpy(T,self.energy,nmax=len(self.evals))

   def getEntropy(self,T):
       return schrodinger.getEntropy(T,self.energy,nmax=len(self.evals))

   def getSumOfStates(self,Elist,sumStates0=None):
       if sumStates0:
           return schrodinger.getSumOfStates(Elist,self.energy,sumStates0=sumStates0,nmax=len(self.evals))
       else:
           return schrodinger.getSumOfStates(Elist,self.energy,nmax=len(self.evals))

   def getDensityOfStates(self,Elist,densStates0=None):
       if densStates0:
           return schrodinger.getDensityOfStates(Elist,self.energy,densStates0=densStates0,nmax=len(self.evals))
       else:
           return schrodinger.getDensityOfStates(Elist,self.energy,nmax=len(self.evals))

   def run(self):
       """
       determines the eigenvalues and energy function for the
       2D-NS rotors either by reading in a finished .evals file
       or running Q2DTor
       """
       try:
           self.readEigvals()
       except IOError:
           self.readScan()
           self.writeXYZ()
           self.getTorsions()
           self.writeInp()
           self.writePes()
           self.getIcsFile()
           self.fitFourier()
           self.getSplistfile()
           self.getEigvals()
