#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import math
import os

import matplotlib.pyplot as plt
import numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import InputError, ElementError, StatmechError
from rmgpy.molecule.molecule import Molecule
from rmgpy.species import TransitionState, Species
from rmgpy.statmech.ndTorsions import HinderedRotor2D, HinderedRotorClassicalND
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor
from rmgpy.statmech.torsion import HinderedRotor, FreeRotor
from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.quantity import Quantity

from arkane.common import ArkaneSpecies, symbol_by_number, get_principal_moments_of_inertia
from arkane.encorr.corr import get_atom_correction, get_bac
from arkane.ess import ESSAdapter, ess_factory, _registered_ess_adapters, GaussianLog, QChemLog
from arkane.encorr.isodesmic import ErrorCancelingSpecies, IsodesmicRingScheme
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory, standardize_name
from arkane.output import prettify
from arkane.encorr.reference import ReferenceDatabase
from arkane.thermo import ThermoJob

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
        'kJ/mol': 1.0 / 1000.,
        'cal/mol': 1.0 / 4.184,
        'kcal/mol': 1.0 / 4184.,
        'cm^-1': 1.0 / (constants.h * constants.c * 100. * constants.Na),
        'hartree': 1.0 / (constants.E_h * constants.Na),
    }

    def __init__(self, path):
        self.path = path

    def load(self):
        """
        Load the scan energies from the file. Returns arrays containing the
        angles (in radians) and energies (in J/mol).
        """
        angles, energies = [], []
        angle_units, energy_units, angle_factor, energy_factor = None, None, None, None

        with open(self.path, 'r') as stream:
            for line in stream:
                line = line.strip()
                if line == '':
                    continue

                tokens = line.split()
                if angle_units is None or energy_units is None:
                    angle_units = tokens[1][1:-1]
                    energy_units = tokens[3][1:-1]

                    try:
                        angle_factor = ScanLog.angleFactors[angle_units]
                    except KeyError:
                        raise ValueError('Invalid angle units {0!r}.'.format(angle_units))
                    try:
                        energy_factor = ScanLog.energyFactors[energy_units]
                    except KeyError:
                        raise ValueError('Invalid energy units {0!r}.'.format(energy_units))

                else:
                    angles.append(float(tokens[0]) / angle_factor)
                    energies.append(float(tokens[1]) / energy_factor)

        angles = np.array(angles)
        energies = np.array(energies)
        energies -= energies[0]

        return angles, energies

    def save(self, angles, energies, angle_units='radians', energy_units='kJ/mol'):
        """
        Save the scan energies to the file using the given `angles` in radians
        and corresponding energies `energies` in J/mol. The file is created to
        use the given `angle_units` for angles and `energy_units` for energies.
        """
        assert len(angles) == len(energies)

        try:
            angle_factor = ScanLog.angleFactors[angle_units]
        except KeyError:
            raise ValueError('Invalid angle units {0!r}.'.format(angle_units))
        try:
            energy_factor = ScanLog.energyFactors[energy_units]
        except KeyError:
            raise ValueError('Invalid energy units {0!r}.'.format(energy_units))

        with open(self.path, 'w') as stream:
            stream.write('{0:>24} {1:>24}\n'.format(
                'Angle ({0})'.format(angle_units),
                'Energy ({0})'.format(energy_units),
            ))
            for angle, energy in zip(angles, energies):
                stream.write('{0:23.10f} {1:23.10f}\n'.format(angle * angle_factor, energy * energy_factor))


################################################################################


def hinderedRotor(scanLog, pivots, top, symmetry=None, fit='best'):
    """Read a hindered rotor directive, and return the attributes in a list"""
    return [scanLog, pivots, top, symmetry, fit]

def hinderedRotor1DArray(angles, energies, pivots, top, symmetry=None, fit='best'):
    """Read a hindered rotor PES profile, and return the attributes in a list"""
    return [angles, energies, pivots, top, symmetry, fit]

def freeRotor(pivots, top, symmetry):
    """Read a free rotor directive, and return the attributes in a list"""
    return [pivots, top, symmetry]


def hinderedRotor2D(scandir, pivots1, top1, symmetry1, pivots2, top2, symmetry2, symmetry='none'):
    """Read a two dimensional hindered rotor directive, and return the attributes in a list"""
    return [scandir, pivots1, top1, symmetry1, pivots2, top2, symmetry2, symmetry]


def hinderedRotorClassicalND(calcPath, pivots, tops, sigmas, semiclassical):
    """Read an N dimensional hindered rotor directive, and return the attributes in a list"""
    return [calcPath, pivots, tops, sigmas, semiclassical]


class StatMechJob(object):
    """
    A representation of a Arkane statistical mechanics job. This job is used
    to compute and save the statistical mechanics information for a single
    species or transition state.
    """

    def __init__(self, species, path):
        self.species = species
        self.path = path
        self.level_of_theory = None
        self.frequencyScaleFactor = 1.0
        self.includeHinderedRotors = True
        self.useIsodesmicReactions = False
        self.isodesmicReactionList = None
        self.referenceSets = None
        self.applyAtomEnergyCorrections = True
        self.applyBondEnergyCorrections = True
        self.bondEnergyCorrectionType = 'p'
        self.atomEnergies = None
        self.bonds = None
        self.arkane_species = ArkaneSpecies(species=species)
        self.hindered_rotor_plots = []

    def execute(self, output_directory=None, plot=False, pdep=False):
        """
        Execute the statmech job, saving the results within
        the `output_directory`.

        If `plot` is True, then plots of the hindered rotor fits will be saved.
        """
        self.load(pdep, plot)
        if output_directory is not None:
            try:
                self.write_output(output_directory)
            except Exception as e:
                logging.warning("Could not write statmech output file due to error: "
                                "{0} in species {1}".format(e, self.species.label))
            if plot:
                hr_dir = os.path.join(output_directory, 'plots')
                if not os.path.exists(hr_dir):
                    os.mkdir(hr_dir)
                try:
                    self.save_hindered_rotor_figures(hr_dir)
                except Exception as e:
                    logging.warning("Could not save hindered rotor scans due to error: "
                                    "{0} in species {1}".format(e, self.species.label))
        logging.debug('Finished statmech job for species {0}.'.format(self.species))
        logging.debug(repr(self.species))

    def load(self, pdep=False, plot=False):
        """
        Load the statistical mechanics parameters for each conformer from
        the associated files on disk. Creates :class:`Conformer` objects for
        each conformer and appends them to the list of conformers on the
        species object.
        """
        path = self.path
        directory = os.path.abspath(os.path.dirname(path))

        def create_log(log_path, check_for_errors=True):
            if not os.path.isfile(log_path):
                modified_log_path = os.path.join(directory, log_path)
                if not os.path.isfile(modified_log_path):
                    raise InputError('Could not find log file for species {0} '
                                     'in the specified path {1}'.format(self.species.label, log_path))
                else:
                    log_path = modified_log_path

            return ess_factory(log_path, check_for_errors=check_for_errors)

        is_ts = isinstance(self.species, TransitionState)
        file_extension = os.path.splitext(path)[-1]
        if file_extension in ['.yml', '.yaml']:
            self.arkane_species.load_yaml(path=path, label=self.species.label, pdep=pdep)
            self.species.conformer = self.arkane_species.conformer
            if is_ts:
                self.species.frequency = self.arkane_species.imaginary_frequency
            else:
                self.species.transport_data = self.arkane_species.transport_data
                self.species.energy_transfer_model = self.arkane_species.energy_transfer_model
                if self.arkane_species.adjacency_list is not None:
                    self.species.molecule = [Molecule().from_adjacency_list(adjlist=self.arkane_species.adjacency_list)]
                elif self.arkane_species.inchi is not None:
                    self.species.molecule = [Molecule().from_inchi(inchistr=self.arkane_species.inchi)]
                elif self.arkane_species.smiles is not None:
                    self.species.molecule = [Molecule().from_smiles(smilesstr=self.arkane_species.smiles)]
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
            'HinderedRotor1DArray': hinderedRotor1DArray,
            'FreeRotor': freeRotor,
            'HinderedRotor2D': hinderedRotor2D,
            'HinderedRotorClassicalND': hinderedRotorClassicalND,
            'ScanLog': ScanLog,
            'Log': create_log,  # The Log class no longer exists, so route the path to ess_factory instead
            'LevelOfTheory': LevelOfTheory,
            'CompositeLevelOfTheory': CompositeLevelOfTheory,
        }

        local_context.update({ess_adapter_name: create_log for ess_adapter_name in _registered_ess_adapters.keys()})

        with open(path, 'r') as f:
            try:
                exec(f.read(), global_context, local_context)
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
            linear = None

        try:
            external_symmetry = local_context['externalSymmetry']
        except KeyError:
            external_symmetry = None

        try:
            spin_multiplicity = local_context['spin_multiplicity']
        except KeyError:
            spin_multiplicity = 0

        try:
            optical_isomers = local_context['opticalIsomers']
        except KeyError:
            logging.debug('No opticalIsomers provided, estimating them from the quantum file.')
            optical_isomers = None

        try:
            energy = local_context['energy']
        except KeyError:
            raise InputError('Required attribute "energy" not found in species file {0!r}.'.format(path))
        if isinstance(energy, dict):
            # Standardize model chemistry names
            energy = {standardize_name(k) if isinstance(k, str) else k: v for k, v in energy.items()}
            freq_level = getattr(self.level_of_theory, 'freq', self.level_of_theory)
            energy_level = getattr(self.level_of_theory, 'energy', self.level_of_theory)
            try:
                energy = energy[energy_level.to_model_chem()]
            except KeyError:
                try:
                    energy = energy[energy_level]
                except KeyError:
                    raise InputError(f'{energy_level} not found in dictionary of energy values in species file {path}.')
        else:
            freq_level = energy_level = None

        e0, e_electronic = None, None  # E0 = e_electronic + ZPE
        energy_log = None
        if isinstance(energy, ESSAdapter):
            energy_log = energy
            # Update energy level of theory with software
            if energy_level is not None:
                energy_software = energy_log.get_software()
                if energy_level.software is not None and energy_level.software != energy_software:
                    logging.warning(f'{energy_level.software} was specified as energy software but does not match'
                                    f' detected software. Software will be updated to {energy_software}.')
                energy_level = energy_level.update(software=energy_software)
        elif isinstance(energy, float):
            e_electronic = energy
        elif isinstance(energy, tuple) and len(energy) == 2:
            # this is likely meant to be a quantity object with ZPE already accounted for
            energy = Quantity(energy)
            e0 = energy.value_si  # in J/mol
        elif isinstance(energy, tuple) and len(energy) == 3:
            if energy[2].lower() == 'e_electronic':
                energy = Quantity(energy[:2])
                e_electronic = energy.value_si / constants.E_h / constants.Na  # convert J/mol to Hartree
            elif energy[2].lower() in ['e0']:
                energy = Quantity(energy[:2])
                e0 = energy.value_si  # in J/mol
            else:
                raise InputError('The third argument for E0 energy value should be e_elect (for energy w/o ZPE) '
                                 'or E0 (including the ZPE). Got: {0}'.format(energy[2]))

        try:
            statmech_log = local_context['frequencies']
        except KeyError:
            raise InputError('Required attribute "frequencies" not found in species file {0!r}.'.format(path))
        try:
            geom_log = local_context['geometry']
        except KeyError:
            geom_log = statmech_log
            logging.debug("Reading geometry from the specified frequencies file.")

        # Update frequency level of theory with software and set new composite level of theory
        if freq_level is not None:
            freq_software = statmech_log.get_software()
            if freq_level.software is not None and freq_level.software != freq_software:
                logging.warning(f'{freq_level.software} was specified as frequency software but does not match detected'
                                f' software. Software will be updated to {freq_software}.')
            freq_level = freq_level.update(software=freq_software)
        if freq_level is not None and energy_level is not None:
            if energy_level == freq_level:
                self.level_of_theory = energy_level
            else:
                self.level_of_theory = CompositeLevelOfTheory(freq=freq_level, energy=energy_level)

        if 'frequencyScaleFactor' in local_context:
            logging.warning('Ignoring frequency scale factor in species file {0!r}.'.format(path))

        rotors = []
        if self.includeHinderedRotors:
            self.raw_hindered_rotor_data = []
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
            if isinstance(statmech_log, GaussianLog):
                if statmech_log.path != geom_log.path:
                    raise InputError('For {0!r}, the geometry log, {1!r}, and frequency log, {2!r}, are not the same. '
                                     'In order to ensure the geometry and Hessian of {0!r} are defined in consistent '
                                     'coordinate systems for hindered/free rotor projection, either use the frequency '
                                     'log for both geometry and frequency, or remove rotors.'.format(
                                      self.species.label, geom_log.path, statmech_log.path))
            elif isinstance(statmech_log, QChemLog):
                logging.warning('QChem log will be used for Hessian of {0!r}. Please verify that the geometry '
                                'and Hessian of {0!r} are defined in the same coordinate system'.format(
                                 self.species.label))

        logging.debug('    Reading molecular degrees of freedom...')
        conformer, unscaled_frequencies = statmech_log.load_conformer(symmetry=external_symmetry,
                                                                      spin_multiplicity=spin_multiplicity,
                                                                      optical_isomers=optical_isomers,
                                                                      label=self.species.label)

        for mode in conformer.modes:
            if isinstance(mode, (Translation, IdealGasTranslation)):
                break
        else:
            # Sometimes the translational mode is not appended to modes for monoatomic species
            conformer.modes.append(IdealGasTranslation(mass=self.species.molecular_weight))

        if conformer.spin_multiplicity == 0:
            raise ValueError("Could not read spin multiplicity from log file {0},\n"
                             "please specify the multiplicity in the input file.".format(self.path))

        logging.debug('    Reading optimized geometry...')
        coordinates, number, mass = geom_log.load_geometry()

        if self.species.conformer is not None and len(self.species.conformer.modes):
            # check that conformer has an IdealGasTranslation mode, append one if it doesn't
            for mode in self.species.conformer.modes:
                if isinstance(mode, IdealGasTranslation):
                    break
            else:
                self.species.conformer.modes.append(IdealGasTranslation(mass=(mass, "amu")))
            # check that conformer has a LinearRotor or a NonlinearRotor mode, append one if it doesn't
            for mode in self.species.conformer.modes:
                if isinstance(mode, (LinearRotor, NonlinearRotor)):
                    break
            else:
                # get the moments of inertia and the external symmetry
                moments_of_inertia = get_principal_moments_of_inertia(coords=self.species.conformer.coordinates,
                                                                      numbers=self.species.conformer.number)
                symmetry = geom_log.get_symmetry_properties()[1]
                if any([moment_of_inertia == 0.0 for moment_of_inertia in moments_of_inertia]):
                    # this is a linear rotor
                    moments_of_inertia = [moment_of_inertia for moment_of_inertia in moments_of_inertia
                                          if moment_of_inertia != 0.0]
                    if abs(moments_of_inertia[0] - moments_of_inertia[1]) > 0.01:
                        raise StatmechError(f'Expected two identical moments of inertia for a linear rigis rotor, '
                                            f'but got {moments_of_inertia}')
                    self.species.conformer.modes.append(LinearRotor(inertia=(moments_of_inertia[0], "amu*angstrom^2"),
                                                                    symmetry=symmetry))
                else:
                    # this is a non-linear rotor
                    self.species.conformer.modes.append(NonlinearRotor(inertia=(moments_of_inertia, "amu*angstrom^2"),
                                                                       symmetry=symmetry))

        # Infer atoms from geometry
        atoms = {}
        for atom_num in number:
            try:
                symbol = symbol_by_number[atom_num]
            except KeyError:
                raise ElementError('Could not recognize element number {0}.'.format(atom_num))
            atoms[symbol] = atoms.get(symbol, 0) + 1

        # Save atoms for use in writing thermo output
        if isinstance(self.species, Species):
            self.species.props['element_counts'] = atoms

        conformer.coordinates = (coordinates, "angstroms")
        conformer.number = number
        conformer.mass = (mass, "amu")

        # The 1.014 factor represents the relationship between the harmonic frequencies scaling factor
        # and the zero point energy scaling factor, see https://pubs.acs.org/doi/10.1021/ct100326h Section 3.1.3.
        zpe_scale_factor = self.frequencyScaleFactor / 1.014

        logging.debug('    Reading energy...')
        if e0 is None:
            if e_electronic is None:
                # The energy read from the log file is without the ZPE
                e_electronic = energy_log.load_energy(zpe_scale_factor)  # in J/mol
            else:
                e_electronic *= constants.E_h * constants.Na  # convert Hartree/particle into J/mol

            # Make sure that isodesmic reactions are configured properly if requested
            if self.useIsodesmicReactions:  # Make sure atom and bond corrections are not applied
                if not self.applyAtomEnergyCorrections:
                    logging.warning('Atom corrections not requested but MUST be used since isodesmic reactions are '
                                    'being used')
                    self.applyAtomEnergyCorrections = True
                if self.applyBondEnergyCorrections:
                    logging.warning('Bond corrections requested but will not be used since isodesmic reactions are '
                                    'being used')
                    self.applyBondEnergyCorrections = False

            # Apply atom corrections
            if self.applyAtomEnergyCorrections:
                atom_corrections = get_atom_correction(self.level_of_theory,
                                                       atoms, self.atomEnergies)

            else:
                atom_corrections = 0
                logging.warning('Atom corrections are not being used. Do not trust energies and thermo.')

            # Apply bond corrections
            if self.applyBondEnergyCorrections:
                if not self.bonds and hasattr(self.species, 'molecule') and self.species.molecule:
                    self.bonds = self.species.molecule[0].enumerate_bonds()
                bond_corrections = get_bac(self.level_of_theory, self.bonds, coordinates, number,
                                           bac_type=self.bondEnergyCorrectionType,
                                           multiplicity=conformer.spin_multiplicity)
            else:
                bond_corrections = 0
            e_electronic_with_corrections = e_electronic + atom_corrections + bond_corrections
            # Get ZPE only for polyatomic species (monoatomic species don't have frequencies, so ZPE = 0)
            zpe = statmech_log.load_zero_point_energy() * zpe_scale_factor if len(number) > 1 else 0
            logging.debug('Scaled zero point energy (ZPE) is {0} J/mol'.format(zpe))

            logging.debug('         Harmonic frequencies scaling factor used = {0:g}'.format(self.frequencyScaleFactor))
            logging.debug('         Zero point energy scaling factor used = {0:g}'.format(zpe_scale_factor))
            logging.debug('         Scaled ZPE (0 K) = {0:g} kcal/mol'.format(zpe / 4184.))

        # If loading a transition state, also read the imaginary frequency
        if is_ts:
            neg_freq = statmech_log.load_negative_frequency()
            self.species.frequency = (neg_freq * self.frequencyScaleFactor, "cm^-1")

        # Read and fit the 1D hindered rotors if applicable
        # If rotors are found, the vibrational frequencies are also
        # recomputed with the torsional modes removed

        hessian = statmech_log.load_force_constant_matrix()

        if hessian is not None and len(mass) > 1 and len(rotors) > 0:
            frequencies, rotors, conformer = self._fit_rotors(rotors, conformer, hessian, is_ts, linear, directory,
                                                              plot)

        elif len(conformer.modes) > 2:
            if len(rotors) > 0:
                logging.warning('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already '
                                'present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for '
                                'Gaussian to generate the force constant matrix, if running Molpro include keyword '
                                'print, hessian')
            frequencies = conformer.modes[2].frequencies.value_si
            rotors = np.array([])
        else:
            if len(rotors) > 0:
                logging.warning('Force Constant Matrix Missing Ignoring rotors, if running Gaussian if not already '
                                'present you need to add the keyword iop(7/33=1) in your Gaussian frequency job for '
                                'Gaussian to generate the force constant matrix, if running Molpro include keyword'
                                'print, hessian')
            frequencies = np.array([])
            rotors = np.array([])

        for mode in conformer.modes:
            if isinstance(mode, HarmonicOscillator):
                mode.frequencies = (frequencies * self.frequencyScaleFactor, "cm^-1")

        if self.useIsodesmicReactions:
            # First, check that a species structure has been given
            if not self.species.molecule:
                raise InputError('A structure must be defined in the species block of the input file to perform '
                                 'isodesmic reaction calculations. For example, append the following to the species '
                                 'block: `structure=SMILES(CC)` using ethane as an example here.')

            # Next, load the reference set database
            reference_db = ReferenceDatabase()
            reference_db.load(paths=self.referenceSets)

            # Set the uncorrected value for E0 on the conformer object so that we can perform the uncorrected thermo job
            conformer.E0 = ((e_electronic_with_corrections + zpe) * 0.001, 'kJ/mol')
            self.species.conformer = conformer

            uncorrected_thermo_job = ThermoJob(species=self.species, thermo_class='nasa')
            uncorrected_thermo_job.generate_thermo()

            uncorrected_thermo = self.species.thermo.get_enthalpy(298)

            scheme = IsodesmicRingScheme(target=ErrorCancelingSpecies(self.species.molecule[0],
                                                                      (uncorrected_thermo, 'J/mol'),
                                                                      self.level_of_theory),
                                         reference_set=reference_db.extract_level_of_theory(self.level_of_theory))
            isodesmic_thermo, self.isodesmicReactionList = scheme.calculate_target_enthalpy()

            # Set the difference as the isodesmic EO correction
            e_electronic_with_corrections += isodesmic_thermo.value_si - uncorrected_thermo

        e0 = e_electronic_with_corrections + zpe
        logging.debug('         E0 (0 K) = {0:g} kcal/mol'.format(e0 / 4184.))
        conformer.E0 = (e0 * 0.001, "kJ/mol")

        # save supporting information for calculation
        self.supporting_info = [self.species.label]
        optical_isomers_read, symmetry_read, point_group_read = statmech_log.get_symmetry_properties()
        self.supporting_info.append(external_symmetry if external_symmetry else symmetry_read)
        self.supporting_info.append(optical_isomers if optical_isomers else optical_isomers_read)
        self.supporting_info.append(point_group_read)
        for mode in conformer.modes:
            if isinstance(mode, (LinearRotor, NonlinearRotor)):
                self.supporting_info.append(mode)
                break
        else:
            self.supporting_info.append(None)
        if unscaled_frequencies:
            self.supporting_info.append(unscaled_frequencies)
        else:
            self.supporting_info.append(None)
        if is_ts:
            self.supporting_info.append(neg_freq)
        else:
            self.supporting_info.append(None)
        self.supporting_info.append(e_electronic)
        self.supporting_info.append(e_electronic + zpe)
        self.supporting_info.append(e0)
        self.supporting_info.append(list([symbol_by_number[x] for x in number]))  # atom symbols
        self.supporting_info.append(coordinates)
        try:
            t1d = energy_log.get_T1_diagnostic()
        except (NotImplementedError, AttributeError):
            t1d = None
        self.supporting_info.append(t1d)
        try:
            d1d = energy_log.get_D1_diagnostic()
        except (NotImplementedError, AttributeError):
            d1d = None
        self.supporting_info.append(d1d)
        # save conformer
        self.species.conformer = conformer

    def _fit_rotors(self, rotors, conformer, hessian, is_ts, linear, directory, plot):
        logging.debug('    Fitting {0} hindered rotors...'.format(len(rotors)))
        rotor_count = 0
        for j, q in enumerate(rotors):
            symmetry = None
            if len(q) == 3:
                # No potential scan is given, this is a free rotor
                pivots, top, symmetry = q
                inertia = conformer.get_internal_reduced_moment_of_inertia(pivots, top) * constants.Na * 1e23
                rotor = FreeRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                conformer.modes.append(rotor)
                rotor_count += 1
            elif len(q) == 8:
                scan_dir, pivots1, top1, symmetry1, pivots2, top2, symmetry2, symmetry = q
                logging.info("Calculating energy levels for 2D-HR, may take a while...")
                rotor = HinderedRotor2D(name='r' + str(j), torsigma1=symmetry1, torsigma2=symmetry2,
                                        symmetry=symmetry, calc_path=os.path.join(directory, scan_dir),
                                        pivots1=pivots1, pivots2=pivots2, top1=top1, top2=top2)
                rotor.run()
                conformer.modes.append(rotor)
                rotor_count += 2
            elif len(q) == 5 and isinstance(q[1][0], list):
                scan_dir, pivots, tops, sigmas, semiclassical = q
                rotor = HinderedRotorClassicalND(pivots, tops, sigmas, calc_path=os.path.join(directory, scan_dir),
                                                 conformer=conformer, F=hessian,
                                                 semiclassical=semiclassical, is_linear=linear, is_ts=is_ts)
                rotor.run()
                conformer.modes.append(rotor)
                rotor_count += len(pivots)
            elif len(q) in [4, 5, 6]:
                # This is a hindered rotor
                if len(q) == 5 and isinstance(q[0], (ESSAdapter, ScanLog)):
                    # A hindered rotor PES from a log file with symmetry assigned
                    scan_log, pivots, top, symmetry, fit = q
                elif len(q) == 5:
                    # A hindered rotor PES from user input arrays with symmetry not assigned
                    # the symmetry number will be derived from the scan
                    angle, v_list, pivots, top, fit = q
                    scan_log = -1
                elif len(q) == 4:
                    # A hindered rotor PES from a log file without symmetry assigned
                    # the symmetry number will be derived from the scan
                    scan_log, pivots, top, fit = q
                elif len(q) == 6:
                    # A hindered rotor PES from user input arrays with symmetry assigned
                    angle, v_list, pivots, top, symmetry, fit = q
                    scan_log = -1
                # Load the hindered rotor scan energies
                if isinstance(scan_log, ScanLog):
                    if not os.path.isfile(scan_log.path):
                        modified_scan_path = os.path.join(directory, scan_log.path)
                        if not os.path.isfile(modified_scan_path):
                            raise InputError('Could not find scan energy file for species {0} '
                                             'in the specified path {1}'.format(self.species.label, scan_log.path))
                        else:
                            scan_log.path = modified_scan_path
                if isinstance(scan_log, (GaussianLog, QChemLog)):
                    v_list, angle = scan_log.load_scan_energies()
                    try:
                        pivot_atoms = scan_log.load_scan_pivot_atoms()
                    except Exception as e:
                        logging.warning("Unable to find pivot atoms in scan due to error: {}".format(e))
                        pivot_atoms = 'N/A'
                    try:
                        frozen_atoms = scan_log.load_scan_frozen_atoms()
                    except Exception as e:
                        logging.warning("Unable to find pivot atoms in scan due to error: {}".format(e))
                        frozen_atoms = 'N/A'
                elif isinstance(scan_log, ScanLog):
                    angle, v_list = scan_log.load()
                    # no way to find pivot atoms or frozen atoms from ScanLog
                    pivot_atoms = 'N/A'
                    frozen_atoms = 'N/A'
                elif scan_log == -1:
                    # Assuming no user may input -1 in the input file. None and '' are not used since they are more likely
                    # to be some input generated from a failure of automatic scripts
                    angle, v_list = np.array(angle), np.array(v_list)
                    pivot_atoms = 'N/A',
                    frozen_atoms = 'N/A'
                else:
                    raise InputError('Invalid log file type {0} for scan log.'.format(scan_log.__class__))

                if symmetry is None:
                    symmetry = determine_rotor_symmetry(v_list, self.species.label, pivots)
                self.raw_hindered_rotor_data.append((self.species.label, rotor_count, symmetry, angle,
                                                     v_list, pivot_atoms, frozen_atoms))
                inertia = conformer.get_internal_reduced_moment_of_inertia(pivots, top) * constants.Na * 1e23

                cosine_rotor = HinderedRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                cosine_rotor.fit_cosine_potential_to_data(angle, v_list)
                fourier_rotor = HinderedRotor(inertia=(inertia, "amu*angstrom^2"), symmetry=symmetry)
                fourier_rotor.fit_fourier_potential_to_data(angle, v_list)

                Vlist_cosine = np.zeros_like(angle)
                Vlist_fourier = np.zeros_like(angle)
                for i in range(angle.shape[0]):
                    Vlist_cosine[i] = cosine_rotor.get_potential(angle[i])
                    Vlist_fourier[i] = fourier_rotor.get_potential(angle[i])

                if fit == 'cosine':
                    rotor = cosine_rotor
                    rotor_count += 1
                    conformer.modes.append(rotor)
                elif fit == 'fourier':
                    rotor = fourier_rotor
                    rotor_count += 1
                    conformer.modes.append(rotor)
                elif fit == 'best':
                    rms_cosine = np.sqrt(np.sum((Vlist_cosine - v_list) * (Vlist_cosine - v_list)) /
                                         (len(v_list) - 1)) / 4184.
                    rms_fourier = np.sqrt(np.sum((Vlist_fourier - v_list) * (Vlist_fourier - v_list)) /
                                          (len(v_list) - 1)) / 4184.

                    # Keep the rotor with the most accurate potential
                    rotor = cosine_rotor if rms_cosine < rms_fourier else fourier_rotor
                    # However, keep the cosine rotor if it is accurate enough, the
                    # fourier rotor is not significantly more accurate, and the cosine
                    # rotor has the correct symmetry
                    if rms_cosine < 0.05 and rms_cosine / rms_fourier < 2.0 and rms_cosine / rms_fourier < 4.0 \
                            and symmetry == cosine_rotor.symmetry:
                        rotor = cosine_rotor

                    conformer.modes.append(rotor)
                    if plot:
                        try:
                            self.create_hindered_rotor_figure(angle, v_list, cosine_rotor, fourier_rotor, rotor,
                                                              rotor_count)
                        except Exception as e:
                            logging.warning("Could not plot hindered rotor graph due to error: {0}".format(e))

                    rotor_count += 1

        logging.debug('    Determining frequencies from reduced force constant matrix...')
        frequencies = np.array(project_rotors(conformer, hessian, rotors, linear, is_ts, label=self.species.label))

        return frequencies, rotors, conformer

    def write_output(self, output_directory):
        """
        Save the results of the statmech job to the `output.py` file located
        in `output_directory`.
        """

        output_file = os.path.join(output_directory, 'output.py')
        logging.info('Saving statistical mechanics parameters for {0}...'.format(self.species.label))
        f = open(output_file, 'a')

        conformer = self.species.conformer
        coordinates = conformer.coordinates.value_si * 1e10
        number = conformer.number.value_si

        f.write('# Coordinates for {0} in Input Orientation (angstroms):\n'.format(self.species.label))
        for i in range(coordinates.shape[0]):
            x = coordinates[i, 0]
            y = coordinates[i, 1]
            z = coordinates[i, 2]
            f.write('#   {0} {1:9.4f} {2:9.4f} {3:9.4f}\n'.format(symbol_by_number[number[i]], x, y, z))

        result = 'conformer(label={0!r}, E0={1!r}, modes={2!r}, spin_multiplicity={3:d}, optical_isomers={4:d}'.format(
            self.species.label,
            conformer.E0,
            conformer.modes,
            conformer.spin_multiplicity,
            conformer.optical_isomers,
        )
        try:
            result += ', frequency={0!r}'.format(self.species.frequency)
        except AttributeError:
            pass
        result += ')'
        f.write('{0}\n\n'.format(prettify(result)))

        if self.useIsodesmicReactions:
            f.write('\n#Isodesmic Reactions Used:\n#------------------------\n#')
            for i, rxn in enumerate(self.isodesmicReactionList):
                thermo = rxn.calculate_target_thermo()
                f.write('Reaction {0}: {1:9.3f} kcal/mol\n#'.format(i+1, thermo.value_si/4184.0))
                reactant_string = '\tReactants:\n#\t\t1.0*{0}\n#'.format(rxn.target.molecule.to_smiles())
                product_string = '\tProducts:\n#'
                for spcs, v in rxn.species.items():
                    if v > 0:  # Product
                        product_string += '\t\t{0}*{1}\n#'.format(v, spcs.molecule.to_smiles())
                    else:  # Reactant
                        reactant_string += '\t\t{0}*{1}\n#'.format(abs(v), spcs.molecule.to_smiles())
                f.write(reactant_string + product_string + '\n#')

            f.write('\n')

        f.close()

    def create_hindered_rotor_figure(self, angle, v_list, cosine_rotor, fourier_rotor, rotor, rotor_index):
        """
        Plot the potential for the rotor, along with its cosine and Fourier
        series potential fits, and save it in the `hindered_rotor_plots` attribute.
        """
        phi = np.arange(0, 6.3, 0.02, np.float64)
        Vlist_cosine = np.zeros_like(phi)
        Vlist_fourier = np.zeros_like(phi)
        for i in range(phi.shape[0]):
            Vlist_cosine[i] = cosine_rotor.get_potential(phi[i])
            Vlist_fourier[i] = fourier_rotor.get_potential(phi[i])

        fig = plt.figure(figsize=(6, 5))
        plt.plot(angle, v_list / 4184., 'ok')
        linespec = '-r' if rotor is cosine_rotor else '--r'
        plt.plot(phi, Vlist_cosine / 4184., linespec)
        linespec = '-b' if rotor is fourier_rotor else '--b'
        plt.plot(phi, Vlist_fourier / 4184., linespec)
        plt.legend(['scan', 'cosine', 'fourier'], loc=1)
        plt.xlim(0, 2 * constants.pi)
        plt.xlabel('Angle')
        plt.ylabel('Potential (kcal/mol)')
        plt.title('{0} hindered rotor #{1:d}'.format(self.species.label, rotor_index + 1))

        axes = fig.get_axes()[0]
        axes.set_xticks([float(j * constants.pi / 4) for j in range(0, 9)])
        axes.set_xticks([float(j * constants.pi / 8) for j in range(0, 17)], minor=True)
        axes.set_xticklabels(
            ['$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'])

        self.hindered_rotor_plots.append((fig, rotor_index))
        plt.close(fig)

    def save_hindered_rotor_figures(self, directory):
        """
        Save hindered rotor plots as set of files of the form
        ``rotor_[species_label]_0.pdf`` in the specified directory
        """
        if hasattr(self, 'hindered_rotor_plots'):
            for fig, rotor_index in self.hindered_rotor_plots:
                fig.savefig(os.path.join(directory, 'rotor_{0}_{1:d}.pdf'.format(self.species.label, rotor_index)))


################################################################################


def is_linear(coordinates):
    """
    Determine whether or not the species is linear from its 3D coordinates
    First, try to reduce the problem into just two dimensions, use 3D if the problem cannot be reduced
    `coordinates` is a numpy.array of the species' xyz coordinates
    """
    # epsilon is in degrees
    # (from our experience, linear molecules have precisely 180.0 degrees between all atom triples)
    epsilon = 0.1

    number_of_atoms = len(coordinates)
    if number_of_atoms == 1:
        return False
    if number_of_atoms == 2:
        return True

    # A tensor containing all distance vectors in the molecule
    d = -np.array([c[:, np.newaxis] - c[np.newaxis, :] for c in coordinates.T])
    for i in range(2, len(coordinates)):
        u1 = d[:, 0, 1] / np.linalg.norm(d[:, 0, 1])  # unit vector between atoms 0 and 1
        u2 = d[:, 1, i] / np.linalg.norm(d[:, 1, i])  # unit vector between atoms 1 and i
        a = math.degrees(np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0)))  # angle between atoms 0, 1, i
        if abs(180 - a) > epsilon and abs(a) > epsilon:
            return False
    return True


def project_rotors(conformer, hessian, rotors, linear, is_ts, get_projected_out_freqs=False, label=None):
    """
    For a given `conformer` with associated force constant matrix `hessian`, lists of
    rotor information `rotors`, `pivots`, and `top1`, and the linearity of the
    molecule `linear`, project out the nonvibrational modes from the force
    constant matrix and use this to determine the vibrational frequencies. The
    list of vibrational frequencies is returned in cm^-1.

    Refer to Gaussian whitepaper (http://gaussian.com/vib/) for procedure to calculate
    harmonic oscillator vibrational frequencies using the force constant matrix.
    """
    n_rotors = 0
    for rotor in rotors:
        if len(rotor) == 8:
            n_rotors += 2
        elif len(rotor) == 5 and isinstance(rotor[1][0], list):
            n_rotors += len(rotor[1])
        else:
            n_rotors += 1

    mass = conformer.mass.value_si
    coordinates = conformer.coordinates.value
    if linear is None:
        linear = is_linear(coordinates)
        if linear:
            logging.info('Determined species {0} to be linear.'.format(label))
    n_atoms = len(conformer.mass.value)
    n_vib = 3 * n_atoms - (5 if linear else 6) - n_rotors - (1 if is_ts else 0)

    # Put origin in center of mass
    xm = 0.0
    ym = 0.0
    zm = 0.0
    totmass = 0.0
    for i in range(n_atoms):
        xm += mass[i] * coordinates[i, 0]
        ym += mass[i] * coordinates[i, 1]
        zm += mass[i] * coordinates[i, 2]
        totmass += mass[i]

    xm /= totmass
    ym /= totmass
    zm /= totmass

    for i in range(n_atoms):
        coordinates[i, 0] -= xm
        coordinates[i, 1] -= ym
        coordinates[i, 2] -= zm
    # Make vector with the root of the mass in amu for each atom
    amass = np.sqrt(mass / constants.amu)

    # Rotation matrix
    inertia = conformer.get_moment_of_inertia_tensor()
    inertia_xyz = np.linalg.eigh(inertia)[1]

    external = 6
    if linear:
        external = 5

    d = np.zeros((n_atoms * 3, external), np.float64)

    # Transform the coordinates to the principal axes
    p = np.dot(coordinates, inertia_xyz)

    for i in range(n_atoms):
        # Projection vectors for translation
        d[3 * i + 0, 0] = amass[i]
        d[3 * i + 1, 1] = amass[i]
        d[3 * i + 2, 2] = amass[i]

    # Construction of the projection vectors for external rotation
    for i in range(n_atoms):
        d[3 * i, 3] = (p[i, 1] * inertia_xyz[0, 2] - p[i, 2] * inertia_xyz[0, 1]) * amass[i]
        d[3 * i + 1, 3] = (p[i, 1] * inertia_xyz[1, 2] - p[i, 2] * inertia_xyz[1, 1]) * amass[i]
        d[3 * i + 2, 3] = (p[i, 1] * inertia_xyz[2, 2] - p[i, 2] * inertia_xyz[2, 1]) * amass[i]
        d[3 * i, 4] = (p[i, 2] * inertia_xyz[0, 0] - p[i, 0] * inertia_xyz[0, 2]) * amass[i]
        d[3 * i + 1, 4] = (p[i, 2] * inertia_xyz[1, 0] - p[i, 0] * inertia_xyz[1, 2]) * amass[i]
        d[3 * i + 2, 4] = (p[i, 2] * inertia_xyz[2, 0] - p[i, 0] * inertia_xyz[2, 2]) * amass[i]
        if not linear:
            d[3 * i, 5] = (p[i, 0] * inertia_xyz[0, 1] - p[i, 1] * inertia_xyz[0, 0]) * amass[i]
            d[3 * i + 1, 5] = (p[i, 0] * inertia_xyz[1, 1] - p[i, 1] * inertia_xyz[1, 0]) * amass[i]
            d[3 * i + 2, 5] = (p[i, 0] * inertia_xyz[2, 1] - p[i, 1] * inertia_xyz[2, 0]) * amass[i]

    # Make sure projection matrix is orthonormal

    inertia = np.identity(n_atoms * 3, np.float64)

    p = np.zeros((n_atoms * 3, 3 * n_atoms + external), np.float64)

    p[:, 0:external] = d[:, 0:external]
    p[:, external:external + 3 * n_atoms] = inertia[:, 0:3 * n_atoms]

    for i in range(3 * n_atoms + external):
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += p[j, i] * p[j, i]
        for j in range(3 * n_atoms):
            if norm > 1E-15:
                p[j, i] /= np.sqrt(norm)
            else:
                p[j, i] = 0.0
        for j in range(i + 1, 3 * n_atoms + external):
            proj = 0.0
            for k in range(3 * n_atoms):
                proj += p[k, i] * p[k, j]
            for k in range(3 * n_atoms):
                p[k, j] -= proj * p[k, i]

    # Order p, there will be vectors that are 0.0
    i = 0
    while i < 3 * n_atoms:
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += p[j, i] * p[j, i]
        if norm < 0.5:
            p[:, i:3 * n_atoms + external - 1] = p[:, i + 1:3 * n_atoms + external]
        else:
            i += 1

    # T is the transformation vector from cartesian to internal coordinates
    T = np.zeros((n_atoms * 3, 3 * n_atoms - external), np.float64)

    T[:, 0:3 * n_atoms - external] = p[:, external:3 * n_atoms]

    # Generate mass-weighted force constant matrix
    # This converts the axes to mass-weighted Cartesian axes
    # Units of Fm are J/m^2*kg = 1/s^2
    weighted_hessian = hessian.copy()
    for i in range(n_atoms):
        for j in range(n_atoms):
            for u in range(3):
                for v in range(3):
                    weighted_hessian[3 * i + u, 3 * j + v] /= math.sqrt(mass[i] * mass[j])

    hessian_int = np.dot(T.T, np.dot(weighted_hessian, T))

    # Get eigenvalues of internal force constant matrix, V = 3N-6 * 3N-6
    eig, v = np.linalg.eigh(hessian_int)

    logging.debug('Frequencies from internal Hessian')
    for i in range(3 * n_atoms - external):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i]) / (2 * math.pi * constants.c * 100))

    # Now we can start thinking about projecting out the internal rotations
    d_int = np.zeros((3 * n_atoms, n_rotors), np.float64)

    counter = 0
    for i, rotor in enumerate(rotors):
        if len(rotor) == 5 and isinstance(rotor[1][0], list):
            scan_dir, pivots_list, tops, sigmas, semiclassical = rotor
        elif len(rotor) == 5 and isinstance(rotor[0], (ESSAdapter, ScanLog)):
            scanLog, pivots, top, symmetry, fit = rotor
            pivots_list, tops = [pivots], [top]
        elif len(rotor) == 5:
            _, _, pivots, top, _ = rotor
            pivots_list, tops = [pivots], [top]
        elif len(rotor) == 3:
            pivots, top, symmetry = rotor
            pivots_list = [pivots]
            tops = [top]
        elif len(rotor) == 8:
            scan_dir, pivots1, top1, symmetry1, pivots2, top2, symmetry2, symmetry = rotor
            pivots_list = [pivots1, pivots2]
            tops = [top1, top2]
        elif len(rotor) == 6:
            _, _, pivots, top, _, _ = rotor
            pivots_list, tops = [pivots], [top]
        else:
            raise ValueError("{} not a proper rotor format".format(rotor))
        for k in range(len(tops)):
            top = tops[k]
            pivots = pivots_list[k]
            # Determine pivot atom
            if pivots[0] in top:
                pivot1 = pivots[0]
                pivot2 = pivots[1]
            elif pivots[1] in top:
                pivot1 = pivots[1]
                pivot2 = pivots[0]
            else:
                raise ValueError('Could not determine pivot atom for rotor {}.'.format(label))
            # Projection vectors for internal rotation
            e12 = coordinates[pivot1 - 1, :] - coordinates[pivot2 - 1, :]
            for j in range(n_atoms):
                atom = j + 1
                if atom in top:
                    e31 = coordinates[atom - 1, :] - coordinates[pivot1 - 1, :]
                    d_int[3 * (atom - 1):3 * (atom - 1) + 3, counter] = np.cross(e31, e12) * amass[atom - 1]
                else:
                    e31 = coordinates[atom - 1, :] - coordinates[pivot2 - 1, :]
                    d_int[3 * (atom - 1):3 * (atom - 1) + 3, counter] = np.cross(e31, -e12) * amass[atom - 1]
            counter += 1

    # Normal modes in mass weighted cartesian coordinates
    vmw = np.dot(T, v)
    eigm = np.zeros((3 * n_atoms - external, 3 * n_atoms - external), np.float64)

    for i in range(3 * n_atoms - external):
        eigm[i, i] = eig[i]

    fm = np.dot(vmw, np.dot(eigm, vmw.T))

    # Internal rotations are not normal modes => project them on the normal modes and orthogonalize
    # d_int_proj =  (3N-6) x (3N) x (3N) x (Nrotors)
    d_int_proj = np.dot(vmw.T, d_int)

    # Reconstruct d_int
    for i in range(n_rotors):
        for j in range(3 * n_atoms):
            d_int[j, i] = 0
            for k in range(3 * n_atoms - external):
                d_int[j, i] += d_int_proj[k, i] * vmw[j, k]

    # Ortho normalize
    for i in range(n_rotors):
        norm = 0.0
        for j in range(3 * n_atoms):
            norm += d_int[j, i] * d_int[j, i]
        for j in range(3 * n_atoms):
            d_int[j, i] /= np.sqrt(norm)
        for j in range(i + 1, n_rotors):
            proj = 0.0
            for k in range(3 * n_atoms):
                proj += d_int[k, i] * d_int[k, j]
            for k in range(3 * n_atoms):
                d_int[k, j] -= proj * d_int[k, i]

    # calculate the frequencies corresponding to the internal rotors
    int_proj = np.dot(fm, d_int)
    kmus = np.array([np.linalg.norm(int_proj[:, i]) for i in range(int_proj.shape[1])])
    int_rotor_freqs = np.sqrt(kmus) / (2.0 * math.pi * constants.c * 100.0)

    if get_projected_out_freqs:
        return int_rotor_freqs

    # Do the projection
    d_int_proj = np.dot(vmw.T, d_int)
    proj = np.dot(d_int, d_int.T)
    inertia = np.identity(n_atoms * 3, np.float64)
    proj = inertia - proj
    fm = np.dot(proj, np.dot(fm, proj))
    # Get eigenvalues of mass-weighted force constant matrix
    eig, v = np.linalg.eigh(fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation

    logging.debug('Frequencies from projected Hessian')
    for i in range(3 * n_atoms):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', r'invalid value encountered in sqrt')
            logging.debug(np.sqrt(eig[i]) / (2 * math.pi * constants.c * 100))

    return np.sqrt(eig[-n_vib:]) / (2 * math.pi * constants.c * 100)


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
        # We declare this rotor as symmetric and the symmetry number is the number of peaks (and valleys)
        symmetry = len(peaks)
        reason = 'number of peaks and valleys, all within the determined resolution criteria'
    if symmetry not in [1, 2, 3]:
        logging.warning('Determined symmetry number {0} for rotor of species {1} between pivots {2}; '
                        'you should make sure this makes sense'.format(symmetry, label, pivots))
    else:
        logging.info('Determined a symmetry number of {0} for rotor of species {1} between pivots {2}'
                     ' based on the {3}.'.format(symmetry, label, pivots, reason))
    return symmetry
