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
This module contains the :class:`ThermoJob` class, used to compute and save the
thermodynamics information for a single species.
"""

import logging
import os.path
import string

import numpy as np

import rmgpy.constants as constants
from rmgpy.chemkin import write_thermo_entry
from rmgpy.exceptions import InputError
from rmgpy.molecule import Molecule
from rmgpy.molecule.util import get_element_count
from rmgpy.species import Species
from rmgpy.statmech.rotation import LinearRotor, NonlinearRotor
from rmgpy.thermo.wilhoit import Wilhoit

from arkane.common import ArkaneSpecies, symbol_by_number
from arkane.output import prettify

################################################################################


class ThermoJob(object):
    """
    A representation of an Arkane thermodynamics job. This job is used to
    compute and save the thermodynamics information for a single species.
    """

    def __init__(self, species, thermo_class):
        self.species = species
        self.thermo_class = thermo_class
        self.arkane_species = ArkaneSpecies(species=species)

    def execute(self, output_directory=None, plot=False):
        """
        Execute the thermodynamics job, saving the results within
        the `output_directory`.

        If `plot` is true, then plots of the raw and fitted values for heat
        capacity, entropy, enthalpy, gibbs free energy, and hindered rotors
        will be saved.
        """
        self.generate_thermo()
        if output_directory is not None:
            try:
                self.write_output(output_directory)
            except Exception as e:
                logging.warning("Could not write output file due to error: "
                                "{0} for species {1}".format(e, self.species.label))
            try:
                self.arkane_species.chemkin_thermo_string = self.write_chemkin(output_directory)
            except Exception as e:
                logging.warning("Could not write chemkin output due to error: "
                                "{0} for species {1}".format(e, self.species.label))
            if self.species.molecule is None or len(self.species.molecule) == 0:
                logging.debug("Not generating a YAML file for species {0}, since its structure wasn't"
                              " specified".format(self.species.label))
            else:
                # We're saving a YAML file for species iff Thermo is called and they're structure is known
                self.arkane_species.update_species_attributes(self.species)
                self.arkane_species.save_yaml(path=output_directory)
            if plot:
                try:
                    self.plot(output_directory)
                except Exception as e:
                    logging.warning("Could not create plots due to error: "
                                    "{0} for species {1}".format(e, self.species.label))

    def generate_thermo(self):
        """
        Generate the thermodynamic data for the species and fit it to the
        desired heat capacity model (as specified in the `thermo_class`
        attribute).
        """
        if self.thermo_class.lower() not in ['wilhoit', 'nasa']:
            raise InputError('Unknown thermodynamic model "{0}".'.format(self.thermo_class))

        species = self.species

        logging.debug('Generating {0} thermo model for {1}...'.format(self.thermo_class, species))

        if species.thermo is not None:
            logging.info("Thermo already generated for species {}. Skipping thermo generation.".format(species))
            return None

        Tlist = np.arange(10.0, 3001.0, 10.0, np.float64)
        Cplist = np.zeros_like(Tlist)
        H298 = 0.0
        S298 = 0.0
        conformer = self.species.conformer
        for i in range(Tlist.shape[0]):
            Cplist[i] += conformer.get_heat_capacity(Tlist[i])
        H298 += conformer.get_enthalpy(298.) + conformer.E0.value_si
        S298 += conformer.get_entropy(298.)

        if not any([isinstance(mode, (LinearRotor, NonlinearRotor)) for mode in conformer.modes]):
            # Monatomic species
            n_freq = 0
            n_rotors = 0
            Cp0 = 2.5 * constants.R
            CpInf = 2.5 * constants.R
        else:
            # Polyatomic species
            linear = True if isinstance(conformer.modes[1], LinearRotor) else False
            n_freq = len(conformer.modes[2].frequencies.value)
            n_rotors = len(conformer.modes[3:])
            Cp0 = (3.5 if linear else 4.0) * constants.R
            CpInf = Cp0 + (n_freq + 0.5 * n_rotors) * constants.R

        wilhoit = Wilhoit()
        if n_freq == 0 and n_rotors == 0:
            wilhoit.Cp0 = (Cplist[0], "J/(mol*K)")
            wilhoit.CpInf = (Cplist[0], "J/(mol*K)")
            wilhoit.B = (500., "K")
            wilhoit.H0 = (0.0, "J/mol")
            wilhoit.S0 = (0.0, "J/(mol*K)")
            wilhoit.H0 = (H298 - wilhoit.get_enthalpy(298.15), "J/mol")
            wilhoit.S0 = (S298 - wilhoit.get_entropy(298.15), "J/(mol*K)")
        else:
            wilhoit.fit_to_data(Tlist, Cplist, Cp0, CpInf, H298, S298, B0=500.0)

        if self.thermo_class.lower() == 'nasa':
            species.thermo = wilhoit.to_nasa(Tmin=10.0, Tmax=3000.0, Tint=500.0)
        else:
            species.thermo = wilhoit

    def write_output(self, output_directory):
        """
        Save the results of the thermodynamics job to the `output.py` file located
        in `output_directory`.
        """
        species = self.species
        output_file = os.path.join(output_directory, 'output.py')
        logging.info('Saving thermo for {0}...'.format(species.label))

        with open(output_file, 'a') as f:
            f.write('# Thermodynamics for {0}:\n'.format(species.label))
            H298 = species.get_thermo_data().get_enthalpy(298) / 4184.
            S298 = species.get_thermo_data().get_entropy(298) / 4.184
            f.write('#   Enthalpy of formation (298 K)   = {0:9.3f} kcal/mol\n'.format(H298))
            f.write('#   Entropy of formation (298 K)    = {0:9.3f} cal/(mol*K)\n'.format(S298))
            f.write('#    =========== =========== =========== =========== ===========\n')
            f.write('#    Temperature Heat cap.   Enthalpy    Entropy     Free energy\n')
            f.write('#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)\n')
            f.write('#    =========== =========== =========== =========== ===========\n')
            for T in [300, 400, 500, 600, 800, 1000, 1500, 2000, 2400]:
                try:
                    Cp = species.get_thermo_data().get_heat_capacity(T) / 4.184
                    H = species.get_thermo_data().get_enthalpy(T) / 4184.
                    S = species.get_thermo_data().get_entropy(T) / 4.184
                    G = species.get_thermo_data().get_free_energy(T) / 4184.
                    f.write('#    {0:11g} {1:11.3f} {2:11.3f} {3:11.3f} {4:11.3f}\n'.format(T, Cp, H, S, G))
                except ValueError:
                    logging.debug("Valid thermo for {0} is outside range for temperature {1}".format(species, T))
            f.write('#    =========== =========== =========== =========== ===========\n')

            thermo_string = 'thermo(label={0!r}, thermo={1!r})'.format(species.label, species.get_thermo_data())
            f.write('{0}\n\n'.format(prettify(thermo_string)))

    def write_chemkin(self, output_directory):
        """
        Appends the thermo block to `chem.inp` and species name to
        `species_dictionary.txt` within the `outut_directory` specified
        """
        species = self.species
        with open(os.path.join(output_directory, 'chem.inp'), 'a') as f:
            if isinstance(species, Species):
                if species.molecule and isinstance(species.molecule[0], Molecule):
                    element_counts = get_element_count(species.molecule[0])
                else:
                    try:
                        element_counts = species.props['element_counts']
                    except KeyError:
                        element_counts = self.element_count_from_conformer()
            else:
                element_counts = {'C': 0, 'H': 0}
            chemkin_thermo_string = write_thermo_entry(species, element_counts=element_counts, verbose=True)
            f.write('{0}\n'.format(chemkin_thermo_string))

        # write species dictionary
        if isinstance(species, Species):
            if species.molecule and isinstance(species.molecule[0], Molecule):
                spec_dict_path = os.path.join(output_directory, 'species_dictionary.txt')
                is_species_in_dict = False
                if os.path.isfile(spec_dict_path):
                    with open(spec_dict_path, 'r') as f:
                        # check whether the species dictionary contains this species, in which case do not re-append
                        for line in f.readlines():
                            if species.label == line.strip():
                                is_species_in_dict = True
                                break
                if not is_species_in_dict:
                    with open(spec_dict_path, 'a') as f:
                        f.write(species.molecule[0].to_adjacency_list(remove_h=False, label=species.label))
                        f.write('\n')
        return chemkin_thermo_string

    def element_count_from_conformer(self):
        """
        Get an element count in a dictionary form (e.g.,  {'C': 3, 'H': 8}) from the species.conformer attribute.

        Returns:
            dict: Element count, keys are element symbols,
                  values are number of occurrences of the element in the molecule.
        """
        element_counts = dict()
        for number in self.species.conformer.number.value_si:
            symbol = symbol_by_number[number]
            if symbol in element_counts:
                element_counts[symbol] += 1
            else:
                element_counts[symbol] = 1
        return element_counts

    def plot(self, output_directory):
        """
        Plot the heat capacity, enthapy, entropy, and Gibbs free energy of the
        fitted thermodynamics model, along with the same values from the
        statistical mechanics model that the thermodynamics model was fitted 
        to. The plot is saved to the file ``thermo.pdf`` in the output
        directory. The plot is not generated if ``matplotlib`` is not installed.
        """
        # Skip this step if matplotlib is not installed
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        Tlist = np.arange(10.0, 2501.0, 10.0)
        Cplist = np.zeros_like(Tlist)
        Cplist1 = np.zeros_like(Tlist)
        Hlist = np.zeros_like(Tlist)
        Hlist1 = np.zeros_like(Tlist)
        Slist = np.zeros_like(Tlist)
        Slist1 = np.zeros_like(Tlist)
        Glist = np.zeros_like(Tlist)
        Glist1 = np.zeros_like(Tlist)

        conformer = self.species.conformer
        thermo = self.species.get_thermo_data()
        for i in range(Tlist.shape[0]):
            try:
                Cplist[i] = conformer.get_heat_capacity(Tlist[i])
                Slist[i] = conformer.get_entropy(Tlist[i])
                Hlist[i] = (conformer.get_enthalpy(Tlist[i]) + conformer.E0.value_si) * 0.001
                Glist[i] = Hlist[i] - Tlist[i] * Slist[i] * 0.001
                Cplist1[i] = thermo.get_heat_capacity(Tlist[i])
                Slist1[i] = thermo.get_entropy(Tlist[i])
                Hlist1[i] = thermo.get_enthalpy(Tlist[i]) * 0.001
                Glist1[i] = thermo.get_free_energy(Tlist[i]) * 0.001
            except (ValueError, AttributeError):
                continue

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle('{0}'.format(self.species.label))
        plt.subplot(2, 2, 1)
        plt.plot(Tlist, Cplist / 4.184, '-r', Tlist, Cplist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Heat capacity (cal/mol*K)')
        plt.legend(['statmech', 'fitted'], loc=4)

        plt.subplot(2, 2, 2)
        plt.plot(Tlist, Slist / 4.184, '-r', Tlist, Slist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Entropy (cal/mol*K)')

        plt.subplot(2, 2, 3)
        plt.plot(Tlist, Hlist / 4.184, '-r', Tlist, Hlist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Enthalpy (kcal/mol)')

        plt.subplot(2, 2, 4)
        plt.plot(Tlist, Glist / 4.184, '-r', Tlist, Glist1 / 4.184, '-b')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Gibbs free energy (kcal/mol)')

        fig.subplots_adjust(left=0.10, bottom=0.08, right=0.95, top=0.95, wspace=0.35, hspace=0.20)

        plot_path = os.path.join(output_directory, 'plots')

        if not os.path.exists(plot_path):
            os.mkdir(plot_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        filename = ''.join(c for c in self.species.label if c in valid_chars) + '.pdf'
        plt.savefig(os.path.join(plot_path, filename))
        plt.close()
