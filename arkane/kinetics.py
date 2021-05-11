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
Arkane kinetics module
"""

import logging
import os.path
import string

import numpy as np

import rmgpy.quantity as quantity
from rmgpy.exceptions import SpeciesError, InputError
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.kinetics.tunneling import Wigner, Eckart
from rmgpy.molecule.draw import MoleculeDrawer, create_new_surface

from arkane.common import ArkaneSpecies
from arkane.output import prettify
from arkane.sensitivity import KineticsSensitivity as SensAnalysis

################################################################################


class KineticsJob(object):
    """
    A representation of an Arkane kinetics job. This job is used to compute
    and save the high-pressure-limit kinetics information for a single reaction.

    `usedTST` - a boolean representing if TST was used to calculate the kinetics
                if kinetics is already given in the input, then it is False.
    `three_params` - a boolean representing if the modified three-parameter Arrhenius equation is used to calculate
                     high pressure kinetic rate coefficients. If it is False, the classical two-parameter Arrhenius
                     equation is used.
    """

    def __init__(self, reaction, Tmin=None, Tmax=None, Tlist=None, Tcount=0, sensitivity_conditions=None,
                 three_params=True):
        self.usedTST = False
        self.Tmin = Tmin if Tmin is not None else (298, 'K')
        self.Tmax = Tmax if Tmax is not None else (2500, 'K')
        self.Tcount = Tcount if Tcount > 3 else 50
        self.three_params = three_params

        if Tlist is not None:
            self.Tlist = Tlist
            self.Tmin = (min(self.Tlist.value_si), 'K')
            self.Tmax = (max(self.Tlist.value_si), 'K')
            self.Tcount = len(self.Tlist.value_si)
        else:
            self.Tlist = (1 / np.linspace(1 / self.Tmax.value_si, 1 / self.Tmin.value_si, self.Tcount), 'K')

        self.reaction = reaction
        self.k_units = None

        if sensitivity_conditions is not None:
            self.sensitivity_conditions = [quantity.Quantity(condition) for condition in sensitivity_conditions]
        else:
            self.sensitivity_conditions = None

        self.arkane_species = ArkaneSpecies(species=self.reaction.transition_state)

    @property
    def Tmin(self):
        """The minimum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmin

    @Tmin.setter
    def Tmin(self, value):
        self._Tmin = quantity.Temperature(value)

    @property
    def Tmax(self):
        """The maximum temperature at which the computed k(T) values are valid, or ``None`` if not defined."""
        return self._Tmax

    @Tmax.setter
    def Tmax(self, value):
        self._Tmax = quantity.Temperature(value)

    @property
    def Tlist(self):
        """The temperatures at which the k(T) values are computed."""
        return self._Tlist

    @Tlist.setter
    def Tlist(self, value):
        self._Tlist = quantity.Temperature(value)

    def execute(self, output_directory=None, plot=True):
        """
        Execute the kinetics job, saving the results within
        the `output_directory`.

        If `plot` is True, then plots of the raw and fitted values for the kinetics
        will be saved.
        """
        self.generate_kinetics()
        if output_directory is not None:
            try:
                self.write_output(output_directory)
            except Exception as e:
                logging.warning("Could not write kinetics output file due to error: "
                                "{0} in reaction {1}".format(e, self.reaction.label))
            try:
                self.write_chemkin(output_directory)
            except Exception as e:
                logging.warning("Could not write kinetics chemkin output due to error: "
                                "{0} in reaction {1}".format(e, self.reaction.label))
            if plot:
                try:
                    self.plot(output_directory)
                except Exception as e:
                    logging.warning("Could not plot kinetics due to error: "
                                    "{0} in reaction {1}".format(e, self.reaction.label))
                try:
                    self.draw(output_directory)
                except Exception as e:
                    logging.warning("Could not draw reaction {1} due to error: {0}".format(e, self.reaction.label))
            if self.sensitivity_conditions is not None:
                logging.info('\n\nRunning sensitivity analysis...')
                SensAnalysis(self, output_directory)
        logging.debug('Finished kinetics job for reaction {0}.'.format(self.reaction))
        logging.debug(repr(self.reaction))

    def generate_kinetics(self):
        """
        Generate the kinetics data for the reaction and fit it to a modified Arrhenius model.
        """

        if isinstance(self.reaction.kinetics, Arrhenius):
            return None
        self.usedTST = True
        kinetics_class = 'Arrhenius'

        tunneling = self.reaction.transition_state.tunneling
        if isinstance(tunneling, Wigner) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transition_state.frequency.value_si, "cm^-1")
        elif isinstance(tunneling, Eckart) and tunneling.frequency is None:
            tunneling.frequency = (self.reaction.transition_state.frequency.value_si, "cm^-1")
            tunneling.E0_reac = (sum([reactant.conformer.E0.value_si
                                      for reactant in self.reaction.reactants]) * 0.001, "kJ/mol")
            tunneling.E0_TS = (self.reaction.transition_state.conformer.E0.value_si * 0.001, "kJ/mol")
            tunneling.E0_prod = (sum([product.conformer.E0.value_si
                                      for product in self.reaction.products]) * 0.001, "kJ/mol")
        elif tunneling is not None:
            if tunneling.frequency is not None:
                # Frequency was given by the user
                pass
            else:
                raise ValueError('Unknown tunneling model {0!r} for reaction {1}.'.format(tunneling, self.reaction))
        logging.debug('Generating {0} kinetics model for {1}...'.format(kinetics_class, self.reaction))
        klist = np.zeros_like(self.Tlist.value_si)
        for i, t in enumerate(self.Tlist.value_si):
            klist[i] = self.reaction.calculate_tst_rate_coefficient(t)
        order = len(self.reaction.reactants)
        klist *= 1e6 ** (order - 1)
        self.k_units = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[order]
        self.K_eq_units = {2: 'mol^2/cm^6', 1: 'mol/cm^3', 0: '       ', -1: 'cm^3/mol', -2: 'cm^6/mol^2'}[
            len(self.reaction.products) - len(self.reaction.reactants)]
        self.k_r_units = {1: 's^-1', 2: 'cm^3/(mol*s)', 3: 'cm^6/(mol^2*s)'}[len(self.reaction.products)]
        self.reaction.kinetics = Arrhenius().fit_to_data(self.Tlist.value_si, klist, kunits=self.k_units,
                                                         three_params=self.three_params)
        self.reaction.elementary_high_p = True

    def write_output(self, output_directory):
        """
        Save the results of the kinetics job to the `output.py` file located
        in `output_directory`.
        """
        reaction = self.reaction

        ks, k0s, k0_revs, k_revs = [], [], [], []

        logging.info('Saving kinetics for {0}...'.format(reaction))

        order = len(self.reaction.reactants)

        factor = 1e6 ** (order - 1)

        f = open(os.path.join(output_directory, 'output.py'), 'a')

        if self.usedTST:
            # If TST is not used, eg. it was given in 'reaction', then this will throw an error.
            f.write('#   ======= =========== =========== =========== ===============\n')
            f.write('#   Temp.   k (TST)     Tunneling   k (TST+T)   Units\n')
            f.write('#   ======= =========== =========== =========== ===============\n')

            if self.Tlist is None:
                t_list = np.array([300, 400, 500, 600, 800, 1000, 1500, 2000])
            else:
                t_list = self.Tlist.value_si

            for T in t_list:
                tunneling = reaction.transition_state.tunneling
                reaction.transition_state.tunneling = None
                try:
                    k0 = reaction.calculate_tst_rate_coefficient(T) * factor
                except SpeciesError:
                    k0 = 0
                reaction.transition_state.tunneling = tunneling
                try:
                    k = reaction.calculate_tst_rate_coefficient(T) * factor
                    kappa = k / k0
                except (SpeciesError, ZeroDivisionError):
                    k = reaction.get_rate_coefficient(T)
                    kappa = 0
                    logging.info("The species in reaction {0} do not have adequate information for TST, "
                                 "using default kinetics values.".format(reaction))
                tunneling = reaction.transition_state.tunneling
                ks.append(k)
                k0s.append(k0)

                f.write('#    {0:4g} K {1:11.3e} {2:11g} {3:11.3e} {4}\n'.format(T, k0, kappa, k, self.k_units))
            f.write('#   ======= =========== =========== =========== ===============\n')
            f.write('\n\n')

            f.write('#   ======= ============ =========== ============ ============= =========\n')
            f.write('#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units\n')
            f.write('#   ======= ============ =========== ============ ============= =========\n')

            # Initialize Object for Converting Units
            if self.K_eq_units != '       ':
                keq_unit_converter = quantity.Units(self.K_eq_units).get_conversion_factor_from_si()
            else:
                keq_unit_converter = 1

            for n, T in enumerate(t_list):
                k = ks[n]
                k0 = k0s[n]
                K_eq = keq_unit_converter * reaction.get_equilibrium_constant(T)  # returns SI units
                k0_rev = k0 / K_eq
                k_rev = k / K_eq
                k0_revs.append(k0_rev)
                k_revs.append(k_rev)
                f.write('#    {0:4g} K {1:11.3e}   {2}  {3:11.3e}   {4:11.3e}      {5}\n'.format(
                    T, K_eq, self.K_eq_units, k0_rev, k_rev, self.k_r_units))

            f.write('#   ======= ============ =========== ============ ============= =========\n')
            f.write('\n\n')

            kinetics_0_rev = Arrhenius().fit_to_data(t_list, np.array(k0_revs), kunits=self.k_r_units,
                                                     three_params=self.three_params)
            kinetics_rev = Arrhenius().fit_to_data(t_list, np.array(k_revs), kunits=self.k_r_units,
                                                   three_params=self.three_params)

            f.write('# k_rev (TST) = {0} \n'.format(kinetics_0_rev))
            f.write('# k_rev (TST+T) = {0} \n\n'.format(kinetics_rev))

        if self.three_params:
            f.write('# kinetics fitted using the modified three-parameter Arrhenius equation '
                    'k = A * (T/T0)^n * exp(-Ea/RT) \n')
        else:
            f.write('# kinetics fitted using the two-parameter Arrhenius equation k = A * exp(-Ea/RT) \n')

        # Reaction path degeneracy is INCLUDED in the kinetics itself!
        rxn_str = 'kinetics(label={0!r}, kinetics={1!r})'.format(reaction.label, reaction.kinetics)
        f.write('{0}\n\n'.format(prettify(rxn_str)))

        f.close()

    def write_chemkin(self, output_directory):
        """
        Appends the kinetics rates to `chem.inp` in `outut_directory`
        """

        # obtain a unit conversion factor
        order = len(self.reaction.reactants)
        factor = 1e6 ** (order - 1)

        reaction = self.reaction
        kinetics = reaction.kinetics
        rxn_str = ''
        if reaction.kinetics.comment:
            for line in reaction.kinetics.comment.split("\n"):
                rxn_str += "! {0}\n".format(line)
        rxn_str += '{0!s:51} {1:9.3e} {2:9.3f} {3:9.3f}\n'.format(
            reaction,
            kinetics.A.value_si * factor,
            kinetics.n.value_si,
            kinetics.Ea.value_si / 4184.,
        )

        with open(os.path.join(output_directory, 'chem.inp'), 'a') as f:
            f.write('{0}\n'.format(rxn_str))

    def save_yaml(self, output_directory):
        """
        Save a YAML file for TSs if structures of the respective reactant/s and product/s are known
        """
        if all([spc.molecule is not None and len(spc.molecule)
                for spc in self.reaction.reactants + self.reaction.products]):
            self.arkane_species.update_species_attributes(self.reaction.transition_state)
            self.arkane_species.reaction_label = self.reaction.label
            self.arkane_species.reactants = [{'label': spc.label, 'adjacency_list': spc.molecule[0].to_adjacency_list()}
                                             for spc in self.reaction.reactants]
            self.arkane_species.products = [{'label': spc.label, 'adjacency_list': spc.molecule[0].to_adjacency_list()}
                                            for spc in self.reaction.products]
            self.arkane_species.save_yaml(path=output_directory)

    def plot(self, output_directory):
        """
        Plot both the raw kinetics data and the Arrhenius fit versus 
        temperature. The plot is saved to the file ``kinetics.pdf`` in the
        output directory. The plot is not generated if ``matplotlib`` is not
        installed.
        """
        import matplotlib.pyplot as plt

        f, ax = plt.subplots()
        if self.Tlist is not None:
            t_list = [t for t in self.Tlist.value_si]
        else:
            t_list = 1000.0 / np.arange(0.4, 3.35, 0.05)
        klist = np.zeros_like(t_list)
        klist2 = np.zeros_like(t_list)
        for i in range(len(t_list)):
            klist[i] = self.reaction.calculate_tst_rate_coefficient(t_list[i])
            klist2[i] = self.reaction.kinetics.get_rate_coefficient(t_list[i])

        order = len(self.reaction.reactants)
        klist *= 1e6 ** (order - 1)
        klist2 *= 1e6 ** (order - 1)
        t_list = [1000.0 / t for t in t_list]
        plt.semilogy(t_list, klist, 'ob', label='TST calculation')
        plt.semilogy(t_list, klist2, '-k', label='Fitted rate')
        plt.legend()
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.reaction.products]))
        plt.title(reaction_str)
        plt.xlabel('1000 / Temperature (K^-1)')
        plt.ylabel('Rate coefficient ({0})'.format(self.k_units))

        plot_path = os.path.join(output_directory, 'plots')

        if not os.path.exists(plot_path):
            os.mkdir(plot_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        plt.savefig(os.path.join(plot_path, filename))
        plt.close()

    def draw(self, output_directory, file_format='pdf'):
        """
        Generate a PDF drawing of the reaction.
        This requires that Cairo and its Python wrapper be available; if not,
        the drawing is not generated.

        You may also generate different formats of drawings, by changing format to
        one of the following: `pdf`, `svg`, `png`.
        """

        drawing_path = os.path.join(output_directory, 'paths')

        if not os.path.exists(drawing_path):
            os.mkdir(drawing_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.reaction.products]))
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        path = os.path.join(drawing_path, filename)

        KineticsDrawer().draw(self.reaction, file_format=file_format, path=path)


class KineticsDrawer(object):
    """
    This class provides functionality for drawing the potential energy surface
    for a high pressure limit reaction using the Cairo 2D graphics engine.
    The most common use case is simply::

        KineticsDrawer().draw(reaction, file_format='png', path='network.png')

    where ``reaction`` is the :class:`Reaction` object to draw. You can also
    pass a dict of options to the constructor to affect how the reaction is drawn.
    """

    def __init__(self, options=None):
        self.options = {
            'structures': True,
            'fontFamily': 'sans',
            'fontSizeNormal': 12,
            'Eunits': 'kJ/mol',
            'padding': 16,
            'wellWidth': 64,
            'wellSpacing': 64,
            'Eslope': 1.5,
            'TSwidth': 16,
            'E0offset': 0.0,
        }
        if options:
            self.options.update(options)
        self.clear()

    def clear(self):
        """Clear the drawer"""
        self.reaction = None
        self.wells = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def _get_energy_range(self):
        """
        Return the minimum and maximum energy in J/mol on the potential energy surface.
        """
        e0_min = min(self.wells[0].E0, self.wells[1].E0, self.reaction.transition_state.conformer.E0.value_si)
        e0_max = max(self.wells[0].E0, self.wells[1].E0, self.reaction.transition_state.conformer.E0.value_si)
        if e0_max - e0_min > 5e5:
            # the energy barrier in one of the reaction directions is larger than 500 kJ/mol, warn the user
            logging.warning('The energy differences between the stationary points of reaction {0} '
                            'seems too large.'.format(self.reaction))
            logging.warning('Got the following energies:\nWell 1: {0} kJ/mol\nTS: {1} kJ/mol\nWell 2: {2}'
                            ' kJ/mol'.format(self.wells[0].E0 / 1000., self.wells[1].E0 / 1000.,
                                             self.reaction.transition_state.conformer.E0.value_si / 1000.))
        return e0_min, e0_max

    def _use_structure_for_label(self, configuration):
        """
        Return ``True`` if the configuration should use molecular structures
        for its labels or ``False`` otherwise.
        """

        # Initialize with the current user option value
        use_structures = self.options['structures']

        # But don't use structures if one or more species in the configuration
        # do not have structure data
        for spec in configuration.species_list:
            if spec.molecule is None or len(spec.molecule) == 0:
                use_structures = False
                break

        return use_structures

    def _get_text_size(self, text, padding=2, file_format='pdf'):
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo

        # Use dummy surface to determine text extents
        surface = create_new_surface(file_format)
        cr = cairo.Context(surface)
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding
        return [0, 0, width, height]

    def _draw_text(self, text, cr, x0, y0, padding=2):
        cr.save()
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        cr.move_to(x0 - extents[0] - padding, y0 - extents[1] + padding)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(text)
        cr.restore()
        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding
        return [0, 0, width, height]

    def _get_label_size(self, configuration, file_format='pdf'):
        width = 0
        height = 0
        bounding_rects = []
        if self._use_structure_for_label(configuration):
            for spec in configuration.species_list:
                rect = MoleculeDrawer().draw(spec.molecule[0], file_format=file_format)[2]
                bounding_rects.append(list(rect))
        else:
            for spec in configuration.species_list:
                bounding_rects.append(self._get_text_size(spec.label, file_format=file_format))

        plus_rect = self._get_text_size('+', file_format=file_format)

        for rect in bounding_rects:
            if width < rect[2]:
                width = rect[2]
            height += rect[3] + plus_rect[3]
        height -= plus_rect[3]

        return [0, 0, width, height]

    def _draw_label(self, configuration, cr, x0, y0, file_format='pdf'):

        bounding_rect = self._get_label_size(configuration, file_format=file_format)
        padding = 2

        use_structures = self._use_structure_for_label(configuration)
        y = y0
        for i, spec in enumerate(configuration.species_list):
            if i > 0:
                rect = self._get_text_size('+', padding=padding, file_format=file_format)
                x = x0 - 0.5 * (rect[2] - bounding_rect[2]) + 2 * padding
                self._draw_text('+', cr, x, y)
                y += rect[3]

            if use_structures:
                molecule_drawer = MoleculeDrawer()
                cr.save()
                rect = molecule_drawer.draw(spec.molecule[0], file_format=file_format)[2]
                cr.restore()
                x = x0 - 0.5 * (rect[2] - bounding_rect[2])
                cr.save()
                molecule_drawer.render(cr, offset=(x, y))
                cr.restore()
                y += rect[3]
            else:
                rect = self._get_text_size(spec.label, padding=padding, file_format=file_format)
                x = x0 - 0.5 * (rect[2] - bounding_rect[2]) + 2 * padding
                self._draw_text(spec.label, cr, x, y)
                y += rect[3]

        return bounding_rect

    def draw(self, reaction, file_format, path=None):
        """
        Draw the potential energy surface for the given `network` as a Cairo
        surface of the given `file_format`. If `path` is given, the surface is
        saved to that location on disk.
        """
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                logging.warning('Cairo not found; potential energy surface will not be drawn.')
                return

        self.reaction = reaction
        self.wells = [Well(self.reaction.reactants), Well(self.reaction.products)]

        # Generate the bounding rectangles for each configuration label
        label_rects = []
        for well in self.wells:
            label_rects.append(self._get_label_size(well, file_format=file_format))

        # Get energy range (use kJ/mol internally)
        e0_min, e0_max = self._get_energy_range()
        e0_min *= 0.001
        e0_max *= 0.001

        # Drawing parameters
        padding = self.options['padding']
        well_width = self.options['wellWidth']
        well_spacing = self.options['wellSpacing']
        e_slope = self.options['Eslope']
        ts_width = self.options['TSwidth']

        e0_offset = self.options['E0offset'] * 0.001

        # Choose multiplier to convert energies to desired units (on figure only)
        e_units = self.options['Eunits']
        try:
            e_mult = {'J/mol': 1.0, 'kJ/mol': 0.001, 'cal/mol': 1.0 / 4.184, 'kcal/mol': 1.0 / 4184.,
                      'cm^-1': 1.0 / 11.962}[e_units]
        except KeyError:
            raise InputError('Invalid value "{0}" for Eunits parameter.'.format(e_units))

        # Determine height required for drawing
        e_height = self._get_text_size('0.0', file_format=file_format)[3] + 6
        y_e0 = (e0_max - 0.0) * e_slope + padding + e_height
        height = (e0_max - e0_min) * e_slope + 2 * padding + e_height + 6
        for i in range(len(self.wells)):
            if 0.001 * self.wells[i].E0 == e0_min:
                height += label_rects[i][3]
                break

        # Determine naive position of each well (one per column)
        coordinates = np.zeros((len(self.wells), 2), np.float64)
        x = padding
        for i in range(len(self.wells)):
            well = self.wells[i]
            rect = label_rects[i]
            this_well_width = max(well_width, rect[2])
            e0 = 0.001 * well.E0
            y = y_e0 - e0 * e_slope
            coordinates[i] = [x + 0.5 * this_well_width, y]
            x += this_well_width + well_spacing
        width = x + padding - well_spacing

        # Determine the rectangles taken up by each well
        # We'll use this to merge columns safely so that wells don't overlap
        well_rects = []
        for i in range(len(self.wells)):
            l, t, w, h = label_rects[i]
            x, y = coordinates[i, :]
            if w < well_width:
                w = well_width
            t -= 6 + e_height
            h += 6 + e_height
            well_rects.append([l + x - 0.5 * w, t + y + 6, w, h])

        # Squish columns together from the left where possible until an isomer is encountered
        old_left = np.min(coordinates[:, 0])
        n_left = - 1
        columns = []
        for i in range(n_left, -1, -1):
            top = well_rects[i][1]
            bottom = top + well_rects[i][3]
            for column in columns:
                for c in column:
                    top0 = well_rects[c][1]
                    bottom0 = top + well_rects[c][3]
                    if (top0 <= top <= bottom0) or (top <= top0 <= bottom):
                        # Can't put it in this column
                        break
                else:
                    # Can put it in this column
                    column.append(i)
                    break
            else:
                # Needs a new column
                columns.append([i])
        for column in columns:
            column_width = max([well_rects[c][2] for c in column])
            x = coordinates[column[0] + 1, 0] - 0.5 * well_rects[column[0] + 1][2] - well_spacing - 0.5 * column_width
            for c in column:
                delta = x - coordinates[c, 0]
                well_rects[c][0] += delta
                coordinates[c, 0] += delta
        new_left = np.min(coordinates[:, 0])
        coordinates[:, 0] -= new_left - old_left

        # Squish columns together from the right where possible until an isomer is encountered
        n_right = 3
        columns = []
        for i in range(n_right, len(self.wells)):
            top = well_rects[i][1]
            bottom = top + well_rects[i][3]
            for column in columns:
                for c in column:
                    top0 = well_rects[c][1]
                    bottom0 = top0 + well_rects[c][3]
                    if (top0 <= top <= bottom0) or (top <= top0 <= bottom):
                        # Can't put it in this column
                        break
                else:
                    # Can put it in this column
                    column.append(i)
                    break
            else:
                # Needs a new column
                columns.append([i])
        for column in columns:
            column_width = max([well_rects[c][2] for c in column])
            x = coordinates[column[0] - 1, 0] + 0.5 * well_rects[column[0] - 1][2] + well_spacing + 0.5 * column_width
            for c in column:
                delta = x - coordinates[c, 0]
                well_rects[c][0] += delta
                coordinates[c, 0] += delta

        width = max([rect[2] + rect[0] for rect in well_rects]) - min([rect[0] for rect in well_rects]) + 2 * padding

        # Draw to the final surface
        surface = create_new_surface(file_format=file_format, target=path, width=width, height=height)
        cr = cairo.Context(surface)

        # Some global settings
        cr.select_font_face("sans")
        cr.set_font_size(self.options['fontSizeNormal'])

        # Fill the background with white
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.paint()
        self._draw_text('E0 ({0})'.format(e_units), cr, 15, 10, padding=2)  # write units

        # Draw reactions
        e0_reac = self.wells[0].E0 * 0.001 - e0_offset
        e0_prod = self.wells[1].E0 * 0.001 - e0_offset
        e0_ts = self.reaction.transition_state.conformer.E0.value_si * 0.001 - e0_offset
        x1, y1 = coordinates[0, :]
        x2, y2 = coordinates[1, :]
        x1 += well_spacing / 2.0
        x2 -= well_spacing / 2.0
        if abs(e0_ts - e0_reac) > 0.1 and abs(e0_ts - e0_prod) > 0.1:
            if len(self.reaction.reactants) == 2:
                if e0_reac < e0_prod:
                    x0 = x1 + well_spacing * 0.5
                else:
                    x0 = x2 - well_spacing * 0.5
            elif len(self.reaction.products) == 2:
                if e0_reac < e0_prod:
                    x0 = x2 - well_spacing * 0.5
                else:
                    x0 = x1 + well_spacing * 0.5
            else:
                x0 = 0.5 * (x1 + x2)
            y0 = y_e0 - (e0_ts + e0_offset) * e_slope
            width1 = (x0 - x1)
            width2 = (x2 - x0)
            # Draw horizontal line for TS
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.set_line_width(2.0)
            cr.move_to(x0 - ts_width / 2.0, y0)
            cr.line_to(x0 + ts_width / 2.0, y0)
            cr.stroke()
            # Add background and text for energy
            e0 = "{0:.1f}".format(e0_ts * 1000. * e_mult)
            extents = cr.text_extents(e0)
            x = x0 - extents[2] / 2.0
            y = y0 - 6.0
            cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            cr.move_to(x, y)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(e0)
            # Draw Bezier curve connecting reactants and products through TS
            cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
            cr.set_line_width(1.0)
            cr.move_to(x1, y1)
            cr.curve_to(x1 + width1 / 8.0, y1, x0 - width1 / 8.0 - ts_width / 2.0, y0, x0 - ts_width / 2.0, y0)
            cr.move_to(x0 + ts_width / 2.0, y0)
            cr.curve_to(x0 + width2 / 8.0 + ts_width / 2.0, y0, x2 - width2 / 8.0, y2, x2, y2)
            cr.stroke()
        else:
            width = (x2 - x1)
            # Draw Bezier curve connecting reactants and products through TS
            cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
            cr.set_line_width(1.0)
            cr.move_to(x1, y1)
            cr.curve_to(x1 + width / 4.0, y1, x2 - width / 4.0, y2, x2, y2)
            cr.stroke()

        # Draw wells (after path reactions so that they are on top)
        for i, well in enumerate(self.wells):
            x0, y0 = coordinates[i, :]
            # Draw horizontal line for well
            cr.set_line_width(4.0)
            cr.move_to(x0 - well_width / 2.0, y0)
            cr.line_to(x0 + well_width / 2.0, y0)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.stroke()
            # Add background and text for energy
            e0 = well.E0 * 0.001 - e0_offset
            e0 = "{0:.1f}".format(e0 * 1000. * e_mult)
            extents = cr.text_extents(e0)
            x = x0 - extents[2] / 2.0
            y = y0 - 6.0
            cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            cr.move_to(x, y)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(e0)
            # Draw background and text for label
            x = x0 - 0.5 * label_rects[i][2]
            y = y0 + 6
            cr.rectangle(x, y, label_rects[i][2], label_rects[i][3])
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            self._draw_label(well, cr, x, y, file_format=file_format)

        # Finish Cairo drawing
        if file_format == 'png':
            surface.write_to_png(path)
        else:
            surface.finish()


class Well(object):
    """
    A helper class representing a "well" of species
    `species_list` is a list of at least one entry
    `E0 `is the sum of all species' E0 in that list
    """

    def __init__(self, species_list):
        self.species_list = species_list
        self.E0 = sum([species.conformer.E0.value_si for species in species_list])
