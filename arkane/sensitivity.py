#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains classes for sensitivity analysis
of kinetics and pressure-dependent jobs.
"""

import logging
import os
import string
import yaml
from copy import deepcopy
import numpy as np

from rmgpy.pdep import Configuration
import rmgpy.quantity as quantity
from rmgpy.species import TransitionState
from rmgpy.kinetics import Arrhenius
from rmgpy.exceptions import InvalidMicrocanonicalRateError, ModifiedStrongCollisionError, PressureDependenceError

################################################################################

# Tolerance (in J/mol) used when checking whether a TransitionState's E0 was synthesized as
# E0(TS) = sum(E0(reactants)) + Ea (as RMG and Arkane's own ILT path do), as opposed to being an
# independent, hand-authored value. 1 kJ/mol is comfortably below the floating-point noise floor
# introduced by unit conversions and energy corrections along that construction path, while still
# being small enough to catch a genuinely independent E0.
SYNTHETIC_TS_E0_TOLERANCE = 1000.0  # J/mol (~1 kJ/mol)


class KineticsSensitivity(object):
    """
    The :class:`KineticsSensitivity` class represents an instance of a sensitivity analysis job
    performed for a KineticsJob. The attributes are:

    =================== ================================================================================================
    Attribute           Description
    =================== ================================================================================================
    `conditions`        A list of the conditions at which the sensitivity coefficients are calculated
    `job`               The KineticsJob object
    `f_rates`           A list of forward rates from `job` at the respective `conditions` in the appropriate units
    `r_rates`           A list of reverse rates from `job` at the respective `conditions` in the appropriate units
    `f_sa_rates`        A dictionary with Species as keys and each value is a list of forward rates from `job` at the
                        respective `conditions` in the appropriate units after perturbing the corresponding Species' E0
    `r_sa_rates`        Same as `f_sa_rates`, only for the reverse direction
    `f_sa_coefficients` A dictionary with Species keys and sensitivity coefficients in the forward direction as values
    `r_sa_coefficients` A dictionary with Species keys and sensitivity coefficients in the reverse direction as values
    =================== ================================================================================================
    """

    def __init__(self, job, output_directory):
        self.job = job
        self.output_directory = output_directory
        self.sensitivity_path = os.path.join(output_directory, 'sensitivity')
        self.conditions = self.job.sensitivity_conditions
        self.f_rates = [self.job.reaction.kinetics.get_rate_coefficient(condition.value_si)
                        for condition in self.conditions]
        kr = self.job.reaction.generate_reverse_rate_coefficient()
        self.r_rates = [kr.get_rate_coefficient(condition.value_si) for condition in self.conditions]
        self.f_sa_rates = {}
        self.r_sa_rates = {}
        self.f_sa_coefficients = {}
        self.r_sa_coefficients = {}
        self.perturbation = quantity.Quantity(1, 'kcal/mol')
        self.execute()

    def execute(self):
        """
        Execute the sensitivity analysis for a :class:KineticsJob: object
        """
        for species in [self.job.reaction.reactants[0], self.job.reaction.products[0],
                        self.job.reaction.transition_state]:
            self.perturb(species)
            self.job.execute(output_file=None, plot=False)  # run the perturbed job
            self.f_sa_rates[species] = [self.job.reaction.kinetics.get_rate_coefficient(condition.value_si)
                                        for condition in self.conditions]
            kr = self.job.reaction.generate_reverse_rate_coefficient()
            self.r_sa_rates[species] = [kr.get_rate_coefficient(condition.value_si)
                                        for condition in self.conditions]
            self.unperturb(species)
            # Calculate the sensitivity coefficients according to dln(r) / dln(E0) = (E0 * dr) / (r * dE0)
            self.f_sa_coefficients[species] = [(self.f_sa_rates[species][i] - self.f_rates[i]) /
                                               (self.perturbation.value_si * self.f_rates[i])
                                               for i in range(len(self.conditions))]
            self.r_sa_coefficients[species] = [(self.r_sa_rates[species][i] - self.r_rates[i]) /
                                               (self.perturbation.value_si * self.r_rates[i])
                                               for i in range(len(self.conditions))]
        self.save()
        self.plot()

    def perturb(self, species):
        """Perturb a species' E0"""
        species.conformer.E0.value_si += self.perturbation.value_si

    def unperturb(self, species):
        """Return the species' E0 to its original value"""
        species.conformer.E0.value_si -= self.perturbation.value_si  # restore E0 to its original value

    def save(self):
        """Save the SA results as tabulated data as well as in YAML format"""
        if not os.path.exists(self.sensitivity_path):
            os.mkdir(self.sensitivity_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reaction_str = '{0} {1} {2}'.format(
            ' + '.join([reactant.label for reactant in self.job.reaction.reactants]),
            '<=>', ' + '.join([product.label for product in self.job.reaction.products]))
        filename = ''.join(c for c in reaction_str if c in valid_chars)
        path = os.path.join(self.sensitivity_path, filename + '.txt')
        sa_data = dict()
        sa_data['structures'] = dict()
        with open(path, 'w') as sa_f:
            sa_f.write("Sensitivity analysis for reaction {0}\n\n"
                       "The semi-normalized sensitivity coefficients are calculated as dln(r)/dE0\n"
                       "by perturbing E0 of each well or TS by {1}, and are given in "
                       "`mol/J` units.\n\n\n".format(reaction_str, self.perturbation))
            for species in self.job.reaction.reactants + self.job.reaction.products:
                if species.label not in sa_data['structures']:
                    sa_data['structures'][species.label] = species.to_adjacency_list()
            reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
            ts_label = self.job.reaction.transition_state.label
            products_label = ' + '.join([reactant.label for reactant in self.job.reaction.products])
            sa_data[reactants_label], sa_data[ts_label], sa_data[products_label] = dict(), dict(), dict()
            max_label = max(len(reactants_label), len(products_label), len(ts_label), 10)
            sa_f.write('========================={0}=============================================\n'
                       '| Direction | Well or TS {1}| Temperature (K) | Sensitivity coefficient |\n'
                       '|-----------+------------{2}+-----------------+-------------------------|\n'
                       .format('=' * (max_label - 10), ' ' * (max_label - 10), '-' * (max_label - 10)))
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    reactants_label, ' ' * (max_label - len(reactants_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]))
                sa_data[reactants_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    products_label, ' ' * (max_label - len(products_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]))
                sa_data[products_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Forward   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    ts_label, ' ' * (max_label - len(ts_label)), condition.value_si,
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]))
                sa_data[ts_label][(condition.value_si, 'K', 'Forward')] = \
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]
            sa_f.write('|-----------+------------{0}+-----------------+-------------------------|\n'.format(
                '-' * (max_label - 10)))
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    reactants_label, ' ' * (max_label - len(reactants_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.reactants[0]][i]))
                sa_data[reactants_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.reactants[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    products_label, ' ' * (max_label - len(products_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.products[0]][i]))
                sa_data[products_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.products[0]][i]
            for i, condition in enumerate(self.conditions):
                sa_f.write('| Reverse   | {0} {1}| {2:6.1f}          | {3:+1.2e}               |\n'.format(
                    ts_label, ' ' * (max_label - len(ts_label)), condition.value_si,
                    self.r_sa_coefficients[self.job.reaction.transition_state][i]))
                sa_data[ts_label][(condition.value_si, 'K', 'Reverse')] = \
                    self.f_sa_coefficients[self.job.reaction.transition_state][i]
            sa_f.write('========================={0}=============================================\n'.format(
                '=' * (max_label - 10)))

        with open(os.path.join(self.sensitivity_path, filename + '.yml'), 'w') as f:
            yaml.dump(data=sa_data, stream=f)

    def plot(self):
        """Plot the SA results as horizontal bars"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
        ts_label = self.job.reaction.transition_state.label
        products_label = ' + '.join([reactant.label for reactant in self.job.reaction.products])

        plt.rcdefaults()
        ax = plt.subplots(nrows=len(self.conditions), ncols=2, tight_layout=True)[1]
        labels = [reactants_label, ts_label, products_label]
        min_sa = min(min(min(self.f_sa_coefficients.values())), min(min(self.r_sa_coefficients.values())))
        max_sa = max(max(max(self.f_sa_coefficients.values())), max(max(self.r_sa_coefficients.values())))
        for i, condition in enumerate(self.conditions):
            f_values = [self.f_sa_coefficients[self.job.reaction.reactants[0]][i],
                        self.f_sa_coefficients[self.job.reaction.transition_state][i],
                        self.f_sa_coefficients[self.job.reaction.products[0]][i]]
            r_values = [self.r_sa_coefficients[self.job.reaction.reactants[0]][i],
                        self.r_sa_coefficients[self.job.reaction.transition_state][i],
                        self.r_sa_coefficients[self.job.reaction.products[0]][i]]
            y_pos = np.arange(3)
            ax[i][0].barh(y_pos, f_values, align='center', color='green')
            ax[i][0].set_yticks(y_pos)
            ax[i][0].set_yticklabels(labels)
            ax[i][0].invert_yaxis()  # labels read top-to-bottom
            ax[i][0].set_xlabel(r'Sensitivity: $\frac{\partial\:\ln{k}}{\partial\:E0}$, ($\frac{J}{mol}$)')
            ax[i][0].set_title('Forward, {0}'.format(condition))
            ax[i][0].set_xlim([min_sa, max_sa])
            ax[i][1].barh(y_pos, r_values, align='center', color='blue')
            ax[i][1].set_yticks(y_pos)
            ax[i][1].set_yticklabels(labels)
            ax[i][1].invert_yaxis()  # labels read top-to-bottom
            ax[i][1].set_xlabel(r'Sensitivity: $\frac{\partial\:\ln{k}}{\partial\:E0}$, ($\frac{J}{mol}$)')
            ax[i][1].set_title('Reverse, {0}'.format(condition))
            ax[i][1].set_xlim([min_sa, max_sa])
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

        if not os.path.exists(self.sensitivity_path):
            os.mkdir(self.sensitivity_path)
        valid_chars = "-_.()<=> %s%s" % (string.ascii_letters, string.digits)
        reactants_label = ' + '.join([reactant.label for reactant in self.job.reaction.reactants])
        products_label = ' + '.join([product.label for product in self.job.reaction.products])
        reaction_str = f'{reactants_label} <=> {products_label}'
        filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
        path = os.path.join(self.sensitivity_path, filename)
        plt.savefig(path)
        plt.close()


def _ts_e0_is_synthetic(rxn, ts_e0_si, ea_si):
    """
    Determine whether a TransitionState's E0 was synthesized from the owning path reaction's
    Arrhenius activation energy as ``E0(TS) = sum(E0(reactants)) + Ea`` (the convention used by
    both RMG, see rmgpy/rmg/pdep.py, and Arkane itself when a TS's E0 is missing, see
    arkane/pdep.py), as opposed to being an independent, physically meaningful value supplied by
    a hand-authored Arkane input (a modeless TS with its own E0).

    `ts_e0_si` and `ea_si` must be the TS's E0 and the reaction's Arrhenius Ea, both in SI units
    (J/mol), taken at the SAME point in a perturb/unperturb pair (i.e. both read before that call's
    own shift is applied to them). Because both quantities are shifted by the identical amount on
    a perturb call and shifted back by the identical amount on the matching unperturb call, their
    difference (ts_e0_si - ea_si) is invariant across the pair, and so is this function's verdict --
    this is essential, since silently disagreeing between the perturb and unperturb calls would
    leave Ea perturbed with no corresponding unperturb (or vice versa), corrupting later iterations.

    A path reaction may be stored in either direction relative to which side the TS sits on, so
    both the reactants-side and the products-side sums are tried before concluding the relation
    does not hold.
    """
    for species_list in (rxn.reactants, rxn.products):
        try:
            e0_sum = sum(spc.conformer.E0.value_si for spc in species_list)
        except AttributeError:
            # A species on this side lacks a conformer/E0 (e.g. a placeholder species without
            # statmech data); treat this side as not matching rather than raising, and fall
            # through to try the other side.
            continue
        if abs((ts_e0_si - ea_si) - e0_sum) <= SYNTHETIC_TS_E0_TOLERANCE:
            return True
    return False


class PDepSensitivity(object):
    """
    The :class:`Sensitivity` class represents an instance of a sensitivity analysis job
    performed for a PressureDependenceJob. The attributes are:

    =================== ================================================================================================
    Attribute           Description
    =================== ================================================================================================
    `conditions`        A list of the conditions (each entry is a list of one T and one P quantities) at which the
                        sensitivity coefficients are calculated
    `job`               The PressureDependenceJob object
    `rates`             A dictionary with net_reactions as keys. Values are lists of forward rates from `job` for the
                        respective path reaction at the respective `conditions` in the appropriate units
    `sa_rates`          A dictionary with string representations of net_reactions as keys. Values are dictionaries with
                        Wells or TransitionStates as keys and each value is a list of forward rates from `job` at the
                        respective `conditions` after perturbing the corresponding well's E0, or, for a TS, its E0
                        and (if the owning path reaction is ILT-based rather than RRKM-based, and the TS's E0 was
                        synthesized from that reaction's Arrhenius Ea in the first place) the Ea of the Arrhenius
                        kinetics used to derive its microcanonical rate, both by the same amount. For such TS rows
                        the reported coefficient is therefore NOT a plain dln(r)/dE0(TS): it is a derivative along
                        the coordinate that raises the synthetic barrier (TS E0 and Ea together)
    `sa_coefficients`   A dictionary with similar structure as `sa_rates`, containing the sensitivity coefficients
                        in the forward direction
    =================== ================================================================================================
    """

    def __init__(self, job, output_directory, perturbation, max_iters=5):
        self.job = job
        self.output_directory = output_directory
        self.sensitivity_path = os.path.join(output_directory, 'sensitivity')
        self.conditions = self.job.sensitivity_conditions
        self.perturbation = perturbation
        
        base_wells = []
        base_wells.extend(self.job.network.reactants)
        base_wells.extend(self.job.network.isomers)
        base_wells.extend(self.job.network.products)
        Emin = min([x.E0 for x in base_wells])
        
        base_transition_states = []
        Emax = Emin
        for rxn in self.job.network.path_reactions:
            base_transition_states.append(rxn.transition_state)
            if rxn.transition_state.conformer.E0.value_si > Emax:
                Emax = rxn.transition_state.conformer.E0.value_si
        if self.perturbation.value_si > 0:
            self.job.network.Emin = Emin
            self.job.network.Emax = Emax + self.perturbation.value_si
        else:
            self.job.network.Emin = Emin + self.perturbation.value_si
            self.job.network.Emax = Emax
        
        self.job.execute(output_file=None, plot=False, print_summary=False)
        
        self.rates = {}
        for rxn in self.job.network.net_reactions:
            self.rates[str(rxn)] = []
            for condition in self.conditions:
                self.rates[str(rxn)].append(rxn.kinetics.get_rate_coefficient(condition[0].value_si,
                                                                              condition[1].value_si))
        self.sa_rates = {}
        self.sa_coefficients = {}
        for rxn in self.job.network.net_reactions:
            self.sa_rates[str(rxn)] = {}
            self.sa_coefficients[str(rxn)] = {}
        
        self.max_iters = max_iters
        self.execute()

    def execute(self):
        """
        Execute the sensitivity analysis for a :class:PressureDependenceJob: object
        """
        base_wells = []
        base_wells.extend(self.job.network.reactants)
        base_wells.extend(self.job.network.isomers)
        base_wells.extend(self.job.network.products)
        base_transition_states = []
        for rxn in self.job.network.path_reactions:
            # if rxn.transition_state is not None:
            base_transition_states.append(rxn.transition_state)
        base_perturbation = deepcopy(self.perturbation)
        base_job = self.job
        for j in range(len(base_wells + base_transition_states)):
            self.perturbation = base_perturbation
            base_entry = (base_wells + base_transition_states)[j]
            failed = False
            c = 0
            wells = []
            transition_states = []
            while c < self.max_iters:
                self.job = deepcopy(base_job)
                wells = []
                wells.extend(self.job.network.reactants)
                wells.extend(self.job.network.isomers)
                wells.extend(self.job.network.products)
                transition_states = []
                for rxn in self.job.network.path_reactions:
                    # if rxn.transition_state is not None:
                    transition_states.append(rxn.transition_state)
                entry = (wells+transition_states)[j]
                if entry in wells:
                    logging.info("\n\nPerturbing well '{0}' by {1}:".format(entry, self.perturbation))
                else:
                    logging.info("\n\nPerturbing TS '{0}' by {1}:".format(entry.label, self.perturbation))
                self.perturb(entry)
                try:
                    self.job.execute(output_file=None, plot=False, print_summary=False)  # run the perturbed job
                    self.unperturb(entry)
                    break
                except (InvalidMicrocanonicalRateError, ModifiedStrongCollisionError) as e:
                    self.unperturb(entry)
                    c += 1
                    self.perturbation = quantity.Quantity(self.perturbation.value/2.0, self.perturbation.units)
                    logging.error("Decreasing perturbation to {}".format(self.perturbation))
                    
            if c == self.max_iters:
                if entry in wells:
                    logging.error("Perturbation of well '{0}' has failed".format(entry))
                else:
                    logging.error("Perturbation of TS '{0}' has failed".format(entry.label))
                failed = True
            for rxn in self.job.network.net_reactions:
                if failed:
                    self.sa_rates[str(rxn)][base_entry] = [np.NaN for condition in self.conditions]
                    self.sa_coefficients[str(rxn)][base_entry] = [np.NaN for i in range(len(self.conditions))]
                else:
                    self.sa_rates[str(rxn)][base_entry] = [rxn.kinetics.get_rate_coefficient(
                        condition[0].value_si, condition[1].value_si) for condition in self.conditions]
                    self.sa_coefficients[str(rxn)][base_entry] = [((self.sa_rates[str(rxn)][base_entry][i]
                                                               - self.rates[str(rxn)][i])) /
                                                             (self.perturbation.value_si * self.rates[str(rxn)][i])
                                                             for i in range(len(self.conditions))]
        self.perturbation = base_perturbation
        self.job = base_job
        self.save(base_wells, base_transition_states)
        self.plot(base_wells, base_transition_states)

    def perturb(self, entry, unperturb=False):
        """
        Perturb E0 of `entry` which could be either a :class:TransitionState or a :class:Configuration
        In the latter case, only the first species in the Configuration.species list is perturbed.
        The perturbation is done by addition of the energy amount in self.perturbation.
        If unperturb is `False`, the perturbation is addition of the energy amount in self.perturbation.
        If unperturb is `False`, this is done by subtracting.

        For a :class:TransitionState belonging to one or more path reactions whose microcanonical
        rate is computed via the inverse Laplace transform (ILT) method rather than RRKM theory
        (i.e. ``not rxn.can_tst()``), perturbing the TS's E0 alone is a structural no-op: the ILT
        convolution and its high-pressure-limit renormalization are driven by the reaction's
        Arrhenius activation energy `Ea`, not by the TS conformer's E0. When, by construction,
        E0(TS) = sum(E0(reactants)) + Ea (see :func:`_ts_e0_is_synthetic`), perturbing `Ea` by the
        same amount as E0 is the physically consistent perturbation, and is applied here in
        addition to the E0 perturbation for every ILT-based path reaction that shares this TS
        object (Arkane input files may point more than one `reaction` block at the same
        `transitionState` label, so a single TS entry can own several path reactions; since the E0
        perturbation above already affects all of them -- the object is shared -- the Ea
        perturbation must match that same blast radius or the two representations of the barrier
        would desynchronize). Hand-authored, modeless transition states whose E0 is an independent
        physical value (i.e. the synthetic relation above does not hold) are left with only their
        E0 perturbed, since moving their Ea as well would not be physically justified. RRKM-based
        (`can_tst()` is `True`) path reactions are unaffected, since RRKM genuinely uses the TS's
        E0 and sum of states, and perturbing both would double-count.
        """
        perturbation = self.perturbation.value_si
        if unperturb:
            perturbation *= -1
        if isinstance(entry, TransitionState):
            # Read the pre-shift TS E0 (and, below, each reaction's pre-shift Ea) so that the
            # synthetic-relation check below is invariant across a perturb/unperturb pair: both
            # E0(TS) and Ea are shifted by the identical `perturbation` here and in the matching
            # call, so their difference -- and hence the check's verdict -- does not change
            # between the two calls. See :func:`_ts_e0_is_synthetic` for the full argument.
            ts_e0_before = entry.conformer.E0.value_si
            entry.conformer.E0 = quantity.Energy(ts_e0_before + perturbation, 'J/mol')
            for rxn in self.job.network.path_reactions:
                if rxn.transition_state is not None and rxn.transition_state is entry and not rxn.can_tst():
                    kinetics = rxn.kinetics if rxn.network_kinetics is None else rxn.network_kinetics
                    if isinstance(kinetics, Arrhenius):
                        ea_before = kinetics.Ea.value_si
                        if _ts_e0_is_synthetic(rxn, ts_e0_before, ea_before):
                            kinetics.Ea = quantity.Energy(ea_before + perturbation, 'J/mol')
                        else:
                            logging.info(
                                "Not perturbing Ea for ILT-based path reaction '{0}': its TS E0 does not "
                                "match sum(E0(reactants)) + Ea (nor the products-side equivalent), so it "
                                "appears to carry an independent, hand-authored E0. Only E0 was "
                                "perturbed.".format(rxn))
                    else:
                        logging.warning(
                            "Not perturbing Ea for ILT-based path reaction '{0}' with kinetics of type "
                            "'{1}': the resulting TS sensitivity coefficient will be meaningless.".format(
                                rxn, kinetics.__class__.__name__))
        elif isinstance(entry, Configuration):
            entry.species[0].conformer.E0 = quantity.Energy(entry.species[0].conformer.E0.value_si + perturbation,
                                                            'J/mol')

    def unperturb(self, entry):
        """A helper function for calling self.perturb cleanly when unperturbing"""
        self.perturb(entry, unperturb=True)

    def save(self, wells, transition_states):
        """Save the SA output as tabulated data as well as in YAML format"""
        if not os.path.exists(os.path.join(self.output_directory, 'sensitivity')):
            os.mkdir(os.path.join(self.output_directory, 'sensitivity'))
        valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
        network_str = self.job.network.label
        filename = os.path.join('sensitivity', ''.join(c for c in network_str if c in valid_chars) + '.txt')
        path = os.path.join(self.output_directory, filename)
        sa_data = dict()
        sa_data['structures'] = dict()
        with open(path, 'w') as sa_f:
            sa_f.write("Sensitivity analysis for network {0}\n\n"
                       "The semi-normalized sensitivity coefficients are calculated as dln(r)/dE0\n"
                       "by perturbing E0 of each well or TS by {1}, and are given in `mol/J` units.\n"
                       "For a TS owned by an ILT-based path reaction whose E0 was synthesized from\n"
                       "that reaction's Arrhenius Ea (E0(TS) = sum(E0(reactants)) + Ea), the Ea is\n"
                       "perturbed by the same amount alongside E0, so that row's coefficient is a\n"
                       "derivative along the coordinate that raises the synthetic barrier (TS E0\n"
                       "and Ea together), not a plain dln(r)/dE0(TS).\n\n\n".format(network_str, self.perturbation))
            for rxn in self.job.network.net_reactions:
                for species in rxn.reactants + rxn.products:
                    if species.label not in sa_data['structures']:
                        sa_data['structures'][species.label] = species.to_adjacency_list()
                reactants_label = ' + '.join([reactant.label for reactant in rxn.reactants])
                products_label = ' + '.join([product.label for product in rxn.products])
                reaction_str = f'{reactants_label} <=> {products_label}'
                sa_data[reaction_str] = dict()
                sa_f.write('  Sensitivity of network reaction ' + reaction_str + ' :' + '\n')
                max_label = 40
                sa_f.write('========================={0}==================================================\n'
                           '| Well or TS {1}| Temperature (K) | Pressure (bar) | Sensitivity coefficient |\n'
                           '|------------{2}+-----------------+----------------+-------------------------|\n'
                           .format('=' * (max_label - 10), ' ' * (max_label - 10), '-' * (max_label - 10)))
                for entry in wells + transition_states:
                    if isinstance(entry, TransitionState):
                        entry_label = '(TS) ' + entry.label
                    elif isinstance(entry, Configuration):
                        entry_label = ' + '.join([species.label for species in entry.species])
                    mod_entry_label = entry_label + ' ' * (max_label - len(entry_label))
                    for i, condition in enumerate(self.conditions):
                        sa_f.write('| {0} | {1:6.1f}          | {2:8.2f}       | {3:+1.2e}               |\n'.format(
                            mod_entry_label, condition[0].value_si, condition[1].value_si * 1e-5,
                            self.sa_coefficients[str(rxn)][entry][i]))
                        condition_tuple = (condition[0].value_si, 'K', condition[1].value_si * 1e-5, 'bar')
                        if condition_tuple not in sa_data[reaction_str]:
                            sa_data[reaction_str][condition_tuple] = dict()
                        sa_data[reaction_str][condition_tuple][entry_label] = self.sa_coefficients[str(rxn)][entry][i]
                sa_f.write('========================={0}=================================================='
                           '\n\n\n'.format('=' * (max_label - 10)))

        with open(os.path.join(self.output_directory, 'sensitivity', 'sa_coefficients.yml'), 'w') as f:
            yaml.dump(data=sa_data, stream=f)

    def plot(self, wells, transition_states):
        """Draw the SA results as horizontal bars"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        for rxn in self.job.network.net_reactions:
            plt.rcdefaults()
            ax = plt.subplots(nrows=len(self.conditions), ncols=1, tight_layout=True)[1]
            labels = [str(entry) for entry in wells]
            labels.extend(ts.label for ts in transition_states)
            max_sa = min_sa = self.sa_coefficients[str(rxn)][wells[0]][0]
            for conformer_sa in self.sa_coefficients[str(rxn)].values():
                for sa_condition in conformer_sa:
                    if min_sa > sa_condition:
                        min_sa = sa_condition
                    if max_sa < sa_condition:
                        max_sa = sa_condition
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
            for i, condition in enumerate(self.conditions):
                values = [self.sa_coefficients[str(rxn)][conf][i] for conf in wells + transition_states]
                y_pos = np.arange(len(labels))
                if len(self.conditions) > 1:
                    axis = ax[i]
                else:
                    axis = ax
                axis.barh(y_pos, values, align='center', color=colors[i % len(colors)])
                axis.set_yticks(y_pos)
                axis.set_yticklabels(labels)
                axis.invert_yaxis()  # labels read top-to-bottom
                # Note: for a TS row belonging to an ILT-based path reaction with a synthetic E0
                # (E0(TS) = sum(E0(reactants)) + Ea), this coefficient is a derivative along the
                # coordinate that raises E0 and Ea together, not a plain dln(k)/dE0(TS); see save().
                axis.set_xlabel(r'Sensitivity: $\frac{\partial\:\ln{k}}{\partial\:E0}$, ($\frac{J}{mol}$)')
                # axis.ticklabel_format('sci')
                axis.set_title('{0}, {1}'.format(condition[0], condition[1]))
                try:
                    axis.set_xlim([min_sa, max_sa])
                except:
                    logging.error("Could not set sensitivity plotting limits may be NaNs.")
                    pass
                axis.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

            if not os.path.exists(self.sensitivity_path):
                os.mkdir(self.sensitivity_path)
            valid_chars = "-_.()<=>+ %s%s" % (string.ascii_letters, string.digits)
            reaction_str = str(rxn)
            filename = ''.join(c for c in reaction_str if c in valid_chars) + '.pdf'
            path = os.path.join(self.sensitivity_path, filename)
            plt.savefig(path)
            plt.close()
